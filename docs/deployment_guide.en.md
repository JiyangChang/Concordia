# Concordia deployment and database preparation guide

English | [中文](deployment_guide.md) (both versions have the same content - keep them in sync when editing)

> For workflow v0.3.22 and later (2026-11-06). Download URLs and software versions were verified on 2026-11-03 / 04.

## 1. Resource planning

| | Whole genome (e.g. a 240 Mb plant genome) | Small test (one chromosome, about 13 Mb) |
|---|---|---|
| CPU | 32–48 cores (a SLURM cluster is better) | 12–14 cores |
| Memory | ≥ 128 GB per node (EDTA, BRAKER4); OMArk about 16 GB | about 16 GB works (with smaller `resources`, and the plant database `Viridiplantae.h5` for OMArk) |
| Database disk | about 130 GB (without FCS-GX); FCS-GX needs another ~470 GB | same (same databases) |
| Result disk | about 100–200 GB (BAMs, EDTA intermediates) | about 10–20 GB |
| conda environments + containers | about 25 GB | same |
| Run time | 2–5 days (cluster, in parallel) | a few hours to 1 day |

Disk advice:
- Put databases and results on a large disk; put conda environments on a local ext4/xfs file system (not NTFS or a
  network drive).
- Databases, results and the container cache of this workflow can live on an NTFS disk; execute permissions and
  symbolic links were verified to work.

## 2. System dependencies

| Software | Purpose | Installation |
|---|---|---|
| Linux x86-64 | — | — |
| conda / mamba (miniforge recommended) | all bioinformatics software | https://github.com/conda-forge/miniforge |
| Apptainer (or Singularity 3.x) | all BRAKER4 tools run in containers | see below |
| git, wget, curl | fetching code and databases | system package manager |

### Apptainer

Ubuntu 22.04 / 24.04 (needs sudo):

```bash
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
apptainer --version
```

Notes:
- **Ubuntu 24.04** restricts unprivileged user namespaces by default (`kernel.apparmor_restrict_unprivileged_userns=1`),
  so Apptainer cannot simply be installed from conda without root; install it from the PPA above, which ships the
  required AppArmor profile.
- **Clusters:** usually `module load apptainer` or `singularity` is available; use it as is.
- **Check:** `apptainer exec docker://alpine:3.20 echo ok` should print `ok`.

## 3. Getting the workflow code

```bash
# BRAKER4 is a git submodule (external/BRAKER4) pinned to the verified commit 6a5ee16; --recursive fetches it too
git clone --recursive <this workflow's repository> concordia && cd annotation
# already cloned but external/BRAKER4 is empty: fetch the submodule
git submodule update --init
```

After updating the workflow code the submodule version may change too; run `git pull && git submodule update` to keep them in step.

## 4. Runtime environment

```bash
# Note: conda/mamba must be on PATH (Snakemake uses it to build the environment of each step), otherwise:
#   "Command 'conda info --json' returned non-zero exit status 127"
export PATH=$HOME/miniforge3/bin:$PATH        # or conda init and log in again
# Environment that runs the workflow itself (Snakemake version as required by BRAKER4)
mamba create -n concordia -c conda-forge -c bioconda \
    python=3.12 snakemake=8.18.2 snakemake-executor-plugin-slurm=2.6.0 pandas pyyaml pytest matplotlib
conda activate concordia

# The software environment of each step is built by Snakemake from workflow/envs/*.yaml (pinned versions).
# Environments go to the conda-prefix directory of profiles/slurm/config.yaml (default ~/annotation_conda_envs,
# about 15-20 GB); when compute nodes have no internet access, build them on the login node first:
cd my_species && snakemake -s /path/to/concordia/workflow/Snakefile -d . --configfile config.yaml \
    --profile /path/to/concordia/profiles/slurm --conda-create-envs-only
```

Settings in `profiles/slurm/config.yaml` to adapt to the cluster:

| Setting | Description |
|---|---|
| `conda-prefix` | Directory of the conda environments. Must be on a shared file system with the same path on the login node and the compute nodes. Environments are named by a hash of the yaml content and of the real path of this directory, so they can be reused by several projects (moving the directory rebuilds all environments; use a symbolic link to change the location). Default `~/annotation_conda_envs` (the home directory is shared on most clusters); if the home directory is too small (about 15–20 GB needed), use a directory on shared storage |
| `default-resources.slurm_partition` / `slurm_account` | partition and accounting account |
| `jobs` | maximum number of jobs queued / running at the same time |

The profiles already enable `software-deployment-method: [conda]`; `--sdm conda` is not needed on the command line.

The configuration template and every profile have an English-commented version (`config/config.example.en.yaml`,
`profiles/*/config.en.yaml`) with exactly the same settings as the Chinese one (`tests/test_config_translations.py`
checks this). Snakemake only reads `config.yaml` in a profile directory; to use the English version, copy
`config.en.yaml` over `config.yaml`.

Environments and versions: `workflow/envs/README.md`. In particular:
- 7 environments are built by default (core, repeat, eviann, agat, qc, omark, function); repeat, eviann, agat and qc have dependency or file conflicts with other software and must stay separate (reasons in `workflow/envs/README.md`).
- InterProScan needs Java 11+, provided by openjdk 11 in the `function` environment, independent of the system Java.

## 5. Databases

### 5.1 One-step download

```bash
bash scripts/download_databases.sh /path/to/databases   # can be re-run; resumes downloads
```

The script handles **one item at a time**: a file is downloaded, verified, decompressed and `sync`ed before the next
one starts; before each step it checks the kernel log and pauses if USB storage reports command timeouts. It is not
parallel because downloading and decompressing several large files at once on a USB disk can make some USB bridge
chips (with the UAS driver) time out repeatedly until the disk stops responding (see section 8). Resulting layout:

| Directory | Content | Download size | Decompressed | Setting |
|---|---|---|---|---|
| `orthodb/Viridiplantae.fa.gz` | OrthoDB v12 plant proteins (BRAKER partition; for other clades use the matching partition, e.g. Metazoa, Fungi) | 3.1 GB | not decompressed | `orthodb_proteins` |
| `orthodb/odb12v2_levels.tab.gz`, `odb12v2_level2species.tab.gz` | OrthoDB level tables: extract a protein subset by level (e.g. Brassicales) to reduce BRAKER4's disk use | 0.3 MB | — | `orthodb_subset.levels_tab`, `orthodb_subset.level2species_tab` |
| `omamer/LUCA.h5` | OMAmer database (for OMArk), needs about 12 GB RAM | 9.9 GB | — | `qc.omark_db` |
| `omamer/Viridiplantae.h5` | plant OMAmer database: replaces LUCA.h5 with < 16 GB RAM (about 2 GB RAM; detects plant contamination only) | 1.5 GB | — | `qc.omark_db` |
| `busco_downloads/lineages/<lineage>` | BUSCO datasets (offline mode); also used by BRAKER4's compleasm | about 0.2 GB each | — | `qc.busco_download_path`, `qc.busco_lineages` |
| `eukaryota_odb12`, `placement_files/`, `file_versions.tsv`, `*.done` in `busco_downloads/lineages/` | extra files needed by BRAKER4's compleasm (`compleasm` step; same layout as compleasm's own downloader) | about 0.2 GB | — | as above |
| `eggnog/` | eggNOG-mapper 5.0.2: `eggnog.db`, `eggnog_proteins.dmnd`, taxa | 12 GB | about 58 GB | `function.eggnog_data_dir` |
| `interproscan/interproscan-5.78-109.0/` | complete InterProScan release (all member databases) | 7.0 GB | about 30–60 GB | `function.interproscan_dir` |
| `dfam/dfam40.0.h5` | Dfam 4.0 FamDB root partition: needed by RepeatMasker 4.2 even with a custom library (EDTA); the workflow links it into the conda environment | 60 MB | 450 MB | `repeat.dfam_dir` |
| `go/go-basic.obo` | GO ontology (GO names and BP/MF/CC categories in the report) | 32 MB | — | `function.go_obo` |
| `uniprot/uniprot_sprot.fasta` | UniProt-SwissProt: used by EviAnn (RNA data only) to decide CDSs (if not set, EviAnn uses the homolog proteins) | 90 MB | about 280 MB | `eviann.swissprot` |
| `kegg/kegg_pathway_names.tsv` | KEGG pathway names (KEGG REST, academic use only) | small | — | `function.kegg_pathway_names` |

### 5.2 Notes on the databases

- **BRAKER4 and the BUSCO datasets:** BRAKER4 runs compleasm as part of gene prediction (CDS hints of BUSCO genes for AUGUSTUS; BUSCO genes dropped by the TSEBRA merge are rescued), which is unrelated to the workflow's BUSCO assessment and cannot be switched off (`skip_busco` only disables BRAKER4's own BUSCO assessment). The workflow writes `qc.busco_download_path` into BRAKER4's `config.ini`, so both share one copy of the datasets, and binds it into the BRAKER4 containers; compleasm refreshes `file_versions.tsv.hash` in that directory on every run, so it must be writable. Without `qc.busco_download_path`, BRAKER4 downloads into the result directory `03_Ab_initio/BRAKER4/busco_downloads` (no longer into `external/BRAKER4`).
- **Choosing BUSCO lineages:** check that the name exists in https://busco-data.ezlab.org/v5/data/file_versions.tsv.
  - Moringa belongs to Brassicales, so use `brassicales_odb12`, plus the broader `embryophyta_odb12`.
  - **odb12 has no `eudicots_odb12`**, only `eudicots_odb10`.
  - The lineage list in the download script can be changed per species.
- **InterProScan:**
  - After unpacking, `python3 setup.py -f interproscan.properties` must be run once to build the HMM indexes. The
    `interproscan_setup` rule does this on the first run and writes a `.pipeline_setup_done` marker into the
    installation directory, so the directory must be writable on the first run (read-only is fine afterwards).
    No separate InterProScan installation is needed: `download_databases.sh` fetches the complete release
    (program + all member databases); point `function.interproscan_dir` to the unpacked directory. Java 11
    comes from the `function` environment.
  - Needs Java 11+.
  - The complete release includes all member databases; `function.interproscan_appl` can restrict the run to some,
    e.g. `Pfam,PANTHER`.
- **eggNOG-mapper:** can also be downloaded with `download_eggnog_data.py --data_dir <dir> -y`; the script fetches the
  same files directly from `http://eggnog5.embl.de/download/emapperdb-5.0.2/`.
- **OrthoDB:** BRAKER4 recommends the partition of the large clade of the target species (Viridiplantae for plants);
  the workflow merges it with the homolog proteins and removes duplicates.
- **KEGG pathway names:** when compute nodes have no internet access, give the path of a downloaded file; with
  `download` the names are fetched at run time.
- **NCBI FCS-GX (optional, off by default):**
  - The database is about 470 GB; at run time it should be in `/dev/shm`, which needs about 500 GB RAM.
  - Download (following the NCBI FCS documentation, after installing bioconda `ncbi-fcs-gx`):
    ```bash
    sync_files.py get --mft https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest --dir /path/to/gxdb
    ```
  - Set `assembly_qc.fcs_gx.db` (e.g. `/path/to/gxdb/all`) and `tax_id` in the configuration.

### 5.3 Online services (not recommended for a whole genome)

| Tool | Online access | Limits |
|---|---|---|
| InterProScan | EBI Job Dispatcher REST API | e-mail address required; at most 30 concurrent jobs; proteins submitted one by one, 30,000 proteins take days |
| eggNOG-mapper | official web server | web submission only, results by e-mail; cannot be embedded in the workflow |

## 6. Configuration

```bash
mkdir my_species && cd my_species           # analysis directory: configuration, sample tables and results; the code directory is only read
cp /path/to/concordia/config/config.example.en.yaml config.yaml
```

Required settings:

| Setting | Description |
|---|---|
| `genome`, `species_prefix` | genome and gene ID prefix |
| `rna.short` / `rna.long` | sample tables (long reads need `platform`: `pacbio_isoseq` / `pacbio_hifi` / `ont_cdna` / `ont_drna`) |
| `homologs` | genomes, GFF3 and proteins of the homologous species |
| `orthodb_proteins` | OrthoDB partition proteins |
| `qc.busco_lineages`, `qc.busco_download_path`, `qc.omark_db` | quality assessment |
| `function.interproscan_dir`, `function.eggnog_data_dir` | functional annotation |
| `braker4.executor` | set to `slurm` on a cluster so that BRAKER4's sub-jobs are submitted to SLURM too |
| `braker4.bind` | paths the containers need (mount points of the databases and the result directory) |

Example: `tests/config.test.yaml` in the repository (mock data, for dry-runs); for organising analysis directories see
"Quick start" in the README.

## 7. Running

```bash
conda activate concordia
cd my_species; A=/path/to/concordia
SM="snakemake -s $A/workflow/Snakefile -d . --configfile config.yaml"
# 1) check
$SM --profile $A/profiles/slurm -n
# 2) run (SLURM; complete example in example_slurm.en.md)
$SM --profile $A/profiles/slurm
# single machine
$SM --profile $A/profiles/local       # environments also in ~/annotation_conda_envs, shared with SLURM
# run up to a step
$SM --profile $A/profiles/local --until evm
$SM --profile $A/profiles/local data_qc
```

- **Resuming:** run the same command again.
- **Debugging a step:** read the `.err` in `<step directory>/logs/`, activate the environment named in the script
  header, and re-run that step alone with `bash <step directory>/scripts/<number>_<rule>.sh`.
- **Complete example:** all steps from setup to submission on a SLURM cluster are in `example_slurm.en.md`.
- **Results:** `<outdir>/05_EVM/final/` holds the final annotation, `<outdir>/07_Report/annotation_report.html` the
  statistics report.

## 8. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| BRAKER4 reports `apptainer: command not found` or cannot create namespaces | Apptainer is not installed, or was installed from conda on Ubuntu 24.04 (blocked by AppArmor). Install it from the PPA as in section 2 |
| GeMoMa prints its tool list and exits | `CLI` was added after the bioconda wrapper; `gemoma.cmd` should be `GeMoMa -Xmx{mem_gb}g` |
| BUSCO cannot find a lineage | The lineage name does not exist (e.g. `eudicots_odb12`), or `busco_download_path` does not point to the directory containing `lineages/` |
| OMArk killed by the system (OOM) | OMAmer loads the whole LUCA.h5 into memory; give `resources.omark.mem_mb` at least 13 GB and make sure it does not run together with other large jobs |
| InterProScan reports a too old Java version | The system Java is 8; the workflow uses openjdk 11 from the `function` environment - check that the profile enables `software-deployment-method: [conda]` (both profiles do by default) and that PATH software is not used instead |
| Building environments fails with `create: The following argument was not expected: --no-default-packages` | Snakemake 8.18.2 is incompatible with mamba 2.x; the profiles set `conda-frontend: conda` (conda uses the libmamba solver by default) |
| BRAKER4 fills the system disk | BRAKER4 copies the protein database about 3 times (run directory, GeneMark-ETP arx/ and data/) and writes ProtHint temporary files: the full OrthoDB plant partition (5.4 GB) needs about 20 GB. With limited disk, use `orthodb_subset.levels` for related levels only (e.g. Brassicales + Malvales, about 2.06 M proteins) |
| OMArk (omamer search) OOM-killed (`Killed` / `Out of memory: Killed process ... (omamer)`) | OMAmer reads the whole database into memory; LUCA.h5 needs about 12 GB. On small machines use `Viridiplantae.h5`, or set `resources.omark.mem_mb` above 14000 so that it does not run together with other jobs |
| AGAT reports `Cannot open agat_tmp_<file name>/result_*.stor` or `Failed to delete agat_tmp_<file name>` | AGAT creates a temporary folder named after the input file in the current directory; two jobs on the same file delete each other's temporary files. The workflow runs every AGAT job in its own temporary directory (`<script>.agat`) |
| InterProScan reports `'sequence' is not an amino acid sequence [...]` | A protein has internal stop codons (written as `.` by gffread). The workflow writes internal stops as X in the functional annotation input and marks these transcripts `internal_stop=yes` in the final GFF3 (`05_EVM/final/internal_stops.tsv`) |
| OMArk hangs at `Downloading taxdump.tar.gz from NCBI FTP site` on its first run, or fails on nodes without internet access | ete3 downloads the NCBI taxonomy on first use (about 80 MB download, 780 MB `~/.etetoolkit/taxa.sqlite`). Run `python -c "from ete3 import NCBITaxa; NCBITaxa()"` once in the omark environment on a node with internet access |
| Error `at least one of rna.short / rna.long must list samples` | No RNA sample was found: check the paths and contents of the `rna.short` / `rna.long` sample tables. Without any RNA-seq, set `rna.allow_none: true` (no-RNA mode, see `rna.allow_none` in the configuration template) |
| AGAT reports `progress bar already finished at .../AGAT/OmniscientI.pm line 390` | AGAT 1.7.0 exits when its progress bar count exceeds the estimated line count (reproducible on GeMoMa output). The workflow writes `00_Data/agat_config.yaml` (`progress_bar: false`) and every AGAT call uses it via `--config` |
| EDTA reports `Error: The RMblast engine is not installed in RepeatMasker!` | Misleading message: the real cause is that RepeatMasker 4.2 cannot find the Dfam FamDB (`FamDB data directory not found`), which the bioconda package does not ship. Set `repeat.dfam_dir`; the `repeat_famdb_setup` step links it into famdb's `Libraries/famdb/` in the environment |
| EDTA reports `Raw TIR results not found`, with TIR-Learner `KeyError: 0` (`check_TIR_TSD.py`) in `.err` | TIR-Learner 3.0.7 is incompatible with pandas 3 (`Series[0]` no longer indexes by position); the `repeat` environment pins `pandas=2.3.3` |
| `calcDivergenceFromAlign.pl` reports `Can't locate RepeatMaskerConfig.pm` | bioconda copies RepeatMasker's utility scripts into `bin/`, and the scripts look for modules relative to their parent directory; the rule points `PERL5LIB` to RepeatMasker's actual installation directory |
| `conda info --json` returns 127 | conda is not on PATH; run `export PATH=$HOME/miniforge3/bin:$PATH` or `conda init` |
| A USB disk freezes under heavy I/O, the kernel log repeatedly shows `uas_eh_abort_handler`, processes hang in state `D` | Some USB bridge chips are unstable with the UAS driver. Fixes: ① use the usb-storage driver: find the enclosure's `<VID>:<PID>` with `lsusb`, `echo "options usb-storage quirks=<VID>:<PID>:u" \| sudo tee /etc/modprobe.d/disable-uas.conf && sudo update-initramfs -u`, effective after a reboot or re-plug; ② plug the disk directly into the computer, not through a dock; ③ keep only databases on the USB disk, data and results on the internal disk; ④ do not decompress several large files in parallel on the USB disk |
| GeneMark licence | GeneMark has licence restrictions for commercial use; check the licence before commercial projects |
