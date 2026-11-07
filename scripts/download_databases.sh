#!/usr/bin/env bash
# Database download for the Concordia genome annotation workflow (URLs verified on 2026-11-04)
#
# Usage:  bash download_databases.sh <target directory>
#
# - One item at a time: a file is downloaded, verified, decompressed and synced before the next one starts
#   (friendly to external / USB disks)
# - Can be re-run: finished items (.done markers in the target directory) are skipped; wget -c resumes downloads
# - Before each step the kernel log is checked: if USB storage reported command timeouts (uas_eh_abort) in the last
#   2 minutes, the script waits for the disk to recover
# - Every download is verified (md5 / gzip CRC / file size); on failure it is deleted and fetched again, up to 3 times
# - Large files are decompressed while computing the md5 of the stream; the written file is then re-read from disk
#   (page cache dropped) and compared, because some USB disks can silently corrupt written data
set -Eeuo pipefail   # -E: the ERR trap also fires inside functions
trap 'echo "[$(date "+%F %T")] ERROR: script exited at line $LINENO (exit code $?)"' ERR
DB=${1:?"usage: bash download_databases.sh <target directory>"}
mkdir -p "$DB" && DB=$(cd "$DB" && pwd) && cd "$DB"

log() { echo "[$(date '+%F %T')] $*"; }
slow() { ionice -c3 nice -n 19 "$@"; }

disk_ok() {   # wait while the USB storage bridge reports aborted commands
    command -v journalctl > /dev/null || return 0
    local n
    while true; do
        n=$(journalctl -k --since "-2 min" --no-pager 2>/dev/null | grep -c uas_eh_abort || true)
        [ "$n" = 0 ] && return 0
        log "USB storage reported $n command timeouts; waiting 60 s for the disk to recover"
        sync; sleep 60
    done
}

fetch() {     # fetch <dir> <url> ...
    local d=$1; shift
    mkdir -p "$d"
    for u in "$@"; do
        disk_ok
        log "downloading $(basename "$u")"
        wget -q -c --tries=5 --waitretry=30 -P "$d" "$u" || { log "ERROR: wget failed: $u"; return 1; }
        sync
    done
}

verified_fetch() {   # verified_fetch <dir> <url> <check: md5:<md5 url>|gzip|size>, up to 3 attempts
    local d=$1 u=$2 check=$3 f i
    f="$d/$(basename "$u")"
    mkdir -p "$d"
    for i in 1 2 3; do
        disk_ok
        log "downloading $(basename "$u") (attempt $i)"
        if ! wget -q -c --tries=5 --waitretry=30 -P "$d" "$u"; then
            log "wget failed (exit code $?); retrying in 30 s"
            sleep 30
            continue
        fi
        sync
        case $check in
            md5:*) if ! wget -q --tries=5 --waitretry=30 -O "$f.md5" "${check#md5:}" || [ ! -s "$f.md5" ]; then
                       log "md5 file download failed; retrying in 30 s (the downloaded data file is kept)"
                       rm -f "$f.md5"; sleep 30; continue
                   fi
                   (cd "$d" && md5sum -c "$(basename "$f").md5") && return 0 ;;
            gzip)  slow gzip -t "$f" && { log "gzip CRC check passed"; return 0; } ;;
            size)  local want have
                   want=$( (curl -sSIL --max-time 60 "$u" || true) | awk 'tolower($1)=="content-length:" {v=$2} END {gsub(/\r/,"",v); print v}')
                   have=$(stat -c %s "$f")
                   [ -n "$want" ] && [ "$want" = "$have" ] && { log "file size matches ($have bytes)"; return 0; } ;;
        esac
        log "verification failed: deleting $(basename "$f") and downloading again"
        rm -f "$f"
    done
    log "ERROR: $(basename "$u") failed verification in all 3 attempts"
    return 1
}

uncached_md5() {   # md5 of a file as stored on disk: drop it from the page cache first, then read it
    ionice -c3 nice -n 19 python3 - "$1" <<'PY'
import hashlib, os, sys
p = sys.argv[1]
fd = os.open(p, os.O_RDONLY)
os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)   # pages are clean after sync -> really evicted
h = hashlib.md5()
while True:
    b = os.read(fd, 8 << 20)
    if not b:
        break
    h.update(b)
os.close(fd)
print(h.hexdigest())
PY
}

write_with_md5() {   # write_with_md5 <out>: stdin -> file; prints the md5 of the stream
    tee "$1" | md5sum | cut -d" " -f1
}

verified_unzip() {   # verified_unzip <file.gz> <out>: decompress and check what was WRITTEN to disk
    local gz=$1 out=$2 want have i
    if [ -s "$out" ]; then          # already decompressed earlier: verify without rewriting
        log "verifying the decompressed $(basename "$out") (md5 of the archive stream vs. md5 of the file on disk)"
        want=$(slow gzip -dc "$gz" | md5sum | cut -d" " -f1)
        sync; have=$(uncached_md5 "$out")
        if [ "$want" = "$have" ]; then
            echo "$want  $(basename "$out")" > "$out.md5"; touch "$out.verified"
            log "$(basename "$out") verified ($have)"; return 0
        fi
        log "$(basename "$out") differs from the archive; decompressing again"
    fi
    for i in 1 2 3; do
        disk_ok
        log "decompressing $(basename "$gz") (attempt $i, computing the md5 of the stream)"
        want=$(slow gzip -dc "$gz" | write_with_md5 "$out.tmp")
        sync
        have=$(uncached_md5 "$out.tmp")
        if [ "$want" = "$have" ]; then
            mv -f "$out.tmp" "$out"
            echo "$want  $(basename "$out")" > "$out.md5"; touch "$out.verified"
            log "$(basename "$out") written and verified ($have)"; return 0
        fi
        log "write verification failed (stream $want, disk $have); deleting and decompressing again"
        rm -f "$out.tmp"
    done
    log "ERROR: $(basename "$out") failed write verification in all 3 attempts"; return 1
}

step() {      # step <name> <function>: run once, mark done
    local name=$1 fn=$2
    if [ -f "$DB/.done.$name" ]; then log "$name already done, skipped"; return; fi
    log "=== $name ==="
    "$fn"
    sync
    touch "$DB/.done.$name"
    log "$name done"
}

dfam() {
    local b=https://www.dfam.org/releases/current/families/FamDB
    fetch dfam "$b/dfam40.0.h5.gz" "$b/dfam40.0.h5.gz.md5"
    (cd dfam && md5sum -c dfam40.0.h5.gz.md5)
    disk_ok; slow gunzip -f dfam/dfam40.0.h5.gz
}

go_kegg() {
    mkdir -p go kegg
    disk_ok; wget -q -c -O go/go-basic.obo https://purl.obolibrary.org/obo/go/go-basic.obo
    curl -sSf --max-time 120 https://rest.kegg.jp/list/pathway -o kegg/kegg_pathway_names.tsv
}

busco() {
    # Adapt the lineage list to the species; names and dates: https://busco-data.ezlab.org/v5/data/file_versions.tsv
    local l
    for l in brassicales_odb12.2026-05-22 embryophyta_odb12.2026-05-22; do
        fetch busco_downloads "https://busco-data.ezlab.org/v5/data/lineages/$l.tar.gz"
        mkdir -p busco_downloads/lineages
        disk_ok; slow tar -xzf "busco_downloads/$l.tar.gz" -C busco_downloads/lineages
        rm -f "busco_downloads/$l.tar.gz"
    done
}

# BRAKER4 runs compleasm for gene prediction (BUSCO CDSpart hints, best_by_compleasm) and shares
# busco_downloads/ with this workflow's BUSCO. compleasm additionally needs, in busco_downloads/lineages/:
# file_versions.tsv(.hash), eukaryota_odb12, the eukaryota placement files, and a <lineage>.done marker for every
# lineage (without a marker it downloads the lineage again). Same files and layout as compleasm's own downloader
# (compleasm_kit/compleasm.py, Downloader), so a later BRAKER4 run downloads nothing.
compleasm() {
    local L=busco_downloads/lineages base=https://busco-data.ezlab.org/v5/data name date md5 cat info f
    mkdir -p "$L/placement_files"
    md5_ok() { [ "$(md5sum < "$1" | cut -d' ' -f1)" = "$2" ] || { log "ERROR: md5 of $(basename "$1") does not match"; return 1; }; }
    rm -f "$L/file_versions.tsv" "$L/file_versions.tsv.hash"      # small; always the current version
    fetch "$L" "$base/file_versions.tsv.hash" "$base/file_versions.tsv"
    md5_ok "$L/file_versions.tsv" "$(awk 'NR==1{print $1}' "$L/file_versions.tsv.hash")"
    touch "$L/file_versions.tsv.done"
    # eukaryota_odb12 (compleasm's default lineage)
    if [ ! -f "$L/eukaryota_odb12.done" ]; then
        read -r name date md5 cat info < <(awk '$1=="eukaryota_odb12" && $5=="lineages"' "$L/file_versions.tsv")
        fetch "$L" "$base/lineages/eukaryota_odb12.$date.tar.gz"
        md5_ok "$L/eukaryota_odb12.$date.tar.gz" "$md5"
        disk_ok; slow tar -xzf "$L/eukaryota_odb12.$date.tar.gz" -C "$L"
        rm -f "$L/eukaryota_odb12.$date.tar.gz"
        touch "$L/eukaryota_odb12.done"
    fi
    # placement files: every eukaryota entry, file name <prefix>.<middle>.<date>.<suffix>.tar.gz
    if [ ! -f "$L/placement_files.done" ]; then
        while read -r name date md5 cat info; do
            f=$(awk -v n="$name" -v d="$date" 'BEGIN{k=split(n,p,"."); m=p[2]; for(i=3;i<k;i++) m=m"."p[i]; print p[1]"."m"."d"."p[k]".tar.gz"}')
            fetch "$L/placement_files" "$base/placement_files/$f"
            md5_ok "$L/placement_files/$f" "$md5"
            slow tar -xzf "$L/placement_files/$f" -C "$L/placement_files"
            rm -f "$L/placement_files/$f"
        done < <(awk '$5=="placement_files" && $1 ~ /eukaryota/' "$L/file_versions.tsv")
        touch "$L/placement_files.done"
    fi
    # markers for the lineages already extracted by the busco step
    for f in "$L"/*_odb*/; do
        [ -f "$f/dataset.cfg" ] && touch "${f%/}.done"
    done
    ls "$L"
}

orthodb() { fetch orthodb https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Viridiplantae.fa.gz; }

orthodb_levels() {   # taxonomy tables for orthodb_subset (BRAKER4 protein subset by OrthoDB level)
    local b=https://data.orthodb.org/v12/download/odb_data_dump f
    for f in odb12v2_levels.tab.gz odb12v2_level2species.tab.gz odb12v2_species.tab.gz; do
        verified_fetch orthodb "$b/$f" gzip
    done
}

omamer() { verified_fetch omamer https://omabrowser.org/All/LUCA.h5 size; }   # no official md5

# Plant OMAmer database (1.5 GB, needs about 2 GB RAM): replaces LUCA.h5 (about 12 GB RAM) on machines with < 16 GB.
# Completeness / consistency are unaffected, but only plant contamination is detected.
omamer_plants() { verified_fetch omamer https://omabrowser.org/All/Viridiplantae.h5 size; }

interproscan() {
    local v=5.78-109.0 b=https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0
    verified_fetch interproscan "$b/interproscan-$v-64-bit.tar.gz" "md5:$b/interproscan-$v-64-bit.tar.gz.md5"
    disk_ok; slow tar -pxzf "interproscan/interproscan-$v-64-bit.tar.gz" -C interproscan
    rm -f "interproscan/interproscan-$v-64-bit.tar.gz"
    # setup.py (HMM indexes) must run once before first use: the workflow's interproscan_setup step does it
}

eggnog() {
    local b=http://eggnog5.embl.de/download/emapperdb-5.0.2 f out
    if [ ! -s eggnog/eggnog.taxa.db ]; then
        verified_fetch eggnog "$b/eggnog.taxa.tar.gz" gzip
        disk_ok; slow tar -xzf eggnog/eggnog.taxa.tar.gz -C eggnog && rm -f eggnog/eggnog.taxa.tar.gz
    fi
    for f in eggnog_proteins.dmnd.gz eggnog.db.gz; do
        out=eggnog/${f%.gz}
        [ -f "$out.verified" ] && { log "$(basename "$out") already verified, skipped"; continue; }
        verified_fetch eggnog "$b/$f" gzip       # gzip CRC of the download
        verified_unzip "eggnog/$f" "$out"        # md5 of the decompressed stream vs. bytes on disk
        rm -f "eggnog/$f"
        sync
    done
}

swissprot() {   # EviAnn protein set (eviann.swissprot) when homolog proteins are not to be used
    local f=uniprot/uniprot_sprot.fasta
    [ -f "$f.verified" ] && { log "uniprot_sprot.fasta already verified, skipped"; return 0; }
    verified_fetch uniprot https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz gzip
    verified_unzip uniprot/uniprot_sprot.fasta.gz "$f"
    rm -f uniprot/uniprot_sprot.fasta.gz
    wget -q -O uniprot/RELEASE.txt https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt || true
}

step dfam dfam
step go_kegg go_kegg
step busco busco
step compleasm compleasm
step orthodb orthodb
step orthodb_levels orthodb_levels
step omamer omamer
step omamer_plants omamer_plants
step interproscan interproscan
step eggnog eggnog
step swissprot swissprot
log "all done"
ls -la "$DB" "$DB"/*/ | sed -n "1,60p" || true   # sed reads everything: no SIGPIPE under pipefail
