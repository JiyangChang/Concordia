#!/usr/bin/env bash
# 基因组注释流程的数据库下载（地址于 2026-11-04 核实）
#
# 用法:  bash download_databases.sh <目标目录> [--limit-rate 15m]
#
# - 逐个下载：一个文件下载、校验、解压、sync 完成后才开始下一个（对 USB 移动硬盘友好）
# - 可重复运行：已完成的项目（目录中有 .done 标记）会跳过；wget -c 断点续传
# - 每步开始前检查内核日志：最近 2 分钟内有 USB 存储命令超时（uas_eh_abort）就等待硬盘恢复
# - 每个文件下载后校验（md5 / gzip CRC / 文件大小），失败则删除重下，最多 3 次
# - 大文件解压限速写入（UNZIP_RATE_MB，默认 15 MB/s，每 256 MB fdatasync），写完后清缓存从盘上重读比对 md5
#   （部分 USB 移动硬盘在 UAS 驱动下可能静默写坏数据，表现为解压后的文件 md5 不符）
set -Eeuo pipefail   # -E: the ERR trap also fires inside functions
trap 'echo "[$(date "+%F %T")] ERROR: 脚本在第 $LINENO 行退出（退出码 $?）"' ERR
DB=${1:?"用法: bash download_databases.sh <目标目录> [--limit-rate 15m]"}
RATE=${3:-15m}
mkdir -p "$DB" && DB=$(cd "$DB" && pwd) && cd "$DB"

log() { echo "[$(date '+%F %T')] $*"; }
slow() { ionice -c3 nice -n 19 "$@"; }

disk_ok() {   # wait while the USB storage bridge reports aborted commands
    command -v journalctl > /dev/null || return 0
    local n
    while true; do
        n=$(journalctl -k --since "-2 min" --no-pager 2>/dev/null | grep -c uas_eh_abort || true)
        [ "$n" = 0 ] && return 0
        log "检测到 $n 次 USB 存储命令超时，暂停 60 秒等待硬盘恢复"
        sync; sleep 60
    done
}

fetch() {     # fetch <dir> <url> ...
    local d=$1; shift
    mkdir -p "$d"
    for u in "$@"; do
        disk_ok
        log "下载 $(basename "$u")"
        wget -q -c --tries=5 --waitretry=30 --limit-rate="$RATE" -P "$d" "$u" || { log "ERROR: wget 失败: $u"; return 1; }
        sync
    done
}

verified_fetch() {   # verified_fetch <dir> <url> <check: md5:<md5 url>|gzip|size>, up to 3 attempts
    local d=$1 u=$2 check=$3 f i
    f="$d/$(basename "$u")"
    mkdir -p "$d"
    for i in 1 2 3; do
        disk_ok
        log "下载 $(basename "$u")（第 $i 次）"
        if ! wget -q -c --tries=5 --waitretry=30 --limit-rate="$RATE" -P "$d" "$u"; then
            log "wget 失败（退出码 $?），30 秒后重试"
            sleep 30
            continue
        fi
        sync
        case $check in
            md5:*) if ! wget -q --tries=5 --waitretry=30 -O "$f.md5" "${check#md5:}" || [ ! -s "$f.md5" ]; then
                       log "md5 文件下载失败，30 秒后重试（不删除已下载的数据文件）"
                       rm -f "$f.md5"; sleep 30; continue
                   fi
                   (cd "$d" && md5sum -c "$(basename "$f").md5") && return 0 ;;
            gzip)  slow gzip -t "$f" && { log "gzip CRC 校验通过"; return 0; } ;;
            size)  local want have
                   want=$( (curl -sSIL --max-time 60 "$u" || true) | awk 'tolower($1)=="content-length:" {v=$2} END {gsub(/\r/,"",v); print v}')
                   have=$(stat -c %s "$f")
                   [ -n "$want" ] && [ "$want" = "$have" ] && { log "文件大小一致（$have 字节）"; return 0; } ;;
        esac
        log "校验失败：删除 $(basename "$f")，重新下载"
        rm -f "$f"
    done
    log "ERROR: $(basename "$u") 三次下载均未通过校验"
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

# python code is passed with -c (not "python3 - <<heredoc"): the heredoc would take over stdin,
# which must carry the decompressed data
THROTTLE_PY=$(cat <<'PY'
import hashlib, os, sys, time
out, rate = sys.argv[1], float(sys.argv[2]) * 1024 * 1024
h, done, since_sync, t0 = hashlib.md5(), 0, 0, time.time()
fd = os.open(out, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
src = sys.stdin.buffer
while True:
    b = src.read(4 << 20)
    if not b:
        break
    h.update(b)
    mv = memoryview(b)
    while mv:
        n = os.write(fd, mv)
        mv = mv[n:]
    done += len(b); since_sync += len(b)
    if since_sync >= 256 << 20:
        os.fdatasync(fd); since_sync = 0
    ahead = done / rate - (time.time() - t0)
    if ahead > 0:
        time.sleep(ahead)
os.fdatasync(fd); os.close(fd)
print(h.hexdigest())
PY
)

throttled_write() {  # throttled_write <out> <rate MB/s>: stdin -> file at a steady rate, fdatasync every 256 MB;
                     # prints the md5 of the stream. Avoids write bursts that made the USB disk drop out.
    ionice -c3 nice -n 19 python3 -c "$THROTTLE_PY" "$1" "$2"
}

verified_unzip() {   # verified_unzip <file.gz> <out>: decompress and check what was WRITTEN to disk
    local gz=$1 out=$2 want have i
    if [ -s "$out" ]; then          # already decompressed earlier: verify without rewriting
        log "校验已解压的 $(basename "$out")（压缩包数据流 md5 对比盘上文件 md5）"
        want=$(slow gzip -dc "$gz" | md5sum | cut -d" " -f1)
        sync; have=$(uncached_md5 "$out")
        if [ "$want" = "$have" ]; then
            echo "$want  $(basename "$out")" > "$out.md5"; touch "$out.verified"
            log "$(basename "$out") 校验通过（$have）"; return 0
        fi
        log "$(basename "$out") 与压缩包不一致，重新解压"
    fi
    for i in 1 2 3; do
        disk_ok
        log "解压 $(basename "$gz")（第 $i 次，限速 ${UNZIP_RATE_MB:-15} MB/s，边写边计算 md5）"
        want=$(slow gzip -dc "$gz" | throttled_write "$out.tmp" "${UNZIP_RATE_MB:-15}")
        sync
        have=$(uncached_md5 "$out.tmp")
        if [ "$want" = "$have" ]; then
            mv -f "$out.tmp" "$out"
            echo "$want  $(basename "$out")" > "$out.md5"; touch "$out.verified"
            log "$(basename "$out") 写入校验通过（$have）"; return 0
        fi
        log "写入校验失败（数据流 $want，盘上 $have），删除后重新解压"
        rm -f "$out.tmp"
    done
    log "ERROR: $(basename "$out") 三次解压均未通过写入校验"; return 1
}

step() {      # step <name> <function>: run once, mark done
    local name=$1 fn=$2
    if [ -f "$DB/.done.$name" ]; then log "$name 已完成，跳过"; return; fi
    log "=== $name ==="
    "$fn"
    sync
    touch "$DB/.done.$name"
    log "$name 完成"
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
    # 按物种修改 lineage 列表；名称和日期见 https://busco-data.ezlab.org/v5/data/file_versions.tsv
    local l
    for l in brassicales_odb12.2026-05-22 embryophyta_odb12.2026-05-22; do
        fetch busco_downloads "https://busco-data.ezlab.org/v5/data/lineages/$l.tar.gz"
        mkdir -p busco_downloads/lineages
        disk_ok; slow tar -xzf "busco_downloads/$l.tar.gz" -C busco_downloads/lineages
        rm -f "busco_downloads/$l.tar.gz"
    done
}

orthodb() { fetch orthodb https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Viridiplantae.fa.gz; }

orthodb_levels() {   # taxonomy tables for orthodb_subset (BRAKER4 protein subset by OrthoDB level)
    local b=https://data.orthodb.org/v12/download/odb_data_dump f
    for f in odb12v2_levels.tab.gz odb12v2_level2species.tab.gz odb12v2_species.tab.gz; do
        verified_fetch orthodb "$b/$f" gzip
    done
}

omamer() { verified_fetch omamer https://omabrowser.org/All/LUCA.h5 size; }   # 官方不提供 md5

# 植物 OMAmer 库（1.5 GB，运行约需 2 GB 内存）：内存不足 16 GB 的机器用它代替 LUCA.h5（约需 12 GB 内存）。
# 完整性/一致性评估不受影响，但只能识别植物来源的污染。
omamer_plants() { verified_fetch omamer https://omabrowser.org/All/Viridiplantae.h5 size; }

interproscan() {
    local v=5.78-109.0 b=https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0
    verified_fetch interproscan "$b/interproscan-$v-64-bit.tar.gz" "md5:$b/interproscan-$v-64-bit.tar.gz.md5"
    disk_ok; slow tar -pxzf "interproscan/interproscan-$v-64-bit.tar.gz" -C interproscan
    rm -f "interproscan/interproscan-$v-64-bit.tar.gz"
    # 首次使用前还需运行 setup.py 建索引：流程中的 interproscan_setup 步骤会自动完成
}

eggnog() {
    local b=http://eggnog5.embl.de/download/emapperdb-5.0.2 f out
    if [ ! -s eggnog/eggnog.taxa.db ]; then
        verified_fetch eggnog "$b/eggnog.taxa.tar.gz" gzip
        disk_ok; slow tar -xzf eggnog/eggnog.taxa.tar.gz -C eggnog && rm -f eggnog/eggnog.taxa.tar.gz
    fi
    for f in eggnog_proteins.dmnd.gz eggnog.db.gz; do
        out=eggnog/${f%.gz}
        [ -f "$out.verified" ] && { log "$(basename "$out") 已校验，跳过"; continue; }
        verified_fetch eggnog "$b/$f" gzip       # gzip CRC of the download
        verified_unzip "eggnog/$f" "$out"        # md5 of the decompressed stream vs. bytes on disk
        rm -f "eggnog/$f"
        sync
    done
}

swissprot() {   # EviAnn protein set (eviann.swissprot) when homolog proteins are not to be used
    local f=uniprot/uniprot_sprot.fasta
    [ -f "$f.verified" ] && { log "uniprot_sprot.fasta 已校验，跳过"; return 0; }
    verified_fetch uniprot https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz gzip
    verified_unzip uniprot/uniprot_sprot.fasta.gz "$f"
    rm -f uniprot/uniprot_sprot.fasta.gz
    wget -q -O uniprot/RELEASE.txt https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt || true
}

step dfam dfam
step go_kegg go_kegg
step busco busco
step orthodb orthodb
step orthodb_levels orthodb_levels
step omamer omamer
step omamer_plants omamer_plants
step interproscan interproscan
step eggnog eggnog
step swissprot swissprot
log "全部完成"
ls -la "$DB" "$DB"/*/ | sed -n "1,60p" || true   # sed reads everything: no SIGPIPE under pipefail
