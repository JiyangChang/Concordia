# Concordia example: running on a SLURM cluster

English | [中文](example_slurm.zh.md) (both versions have the same content - keep them in sync when editing)

This example assumes: compute nodes with 48 cores / 192 GB RAM; the login node has internet access, compute nodes
may not; `/shared` is storage visible from all nodes; `myname` stands for the user name. Replace paths, partition and
account with those of your cluster. Deployment details (database contents, troubleshooting) are in `deployment_guide.md`.

## 1. Where the results go

**Results are written to Snakemake's working directory, not to the directory of the workflow code.** The working
directory is set with `-d`; without `-d` it is the directory from which `snakemake` is run. Relative paths in the
configuration (including `outdir`, sample tables, genome) are resolved against the working directory as well; `~` is
expanded to the home directory, but environment variables such as `$USER` are not - write the actual path.

In this example the working directory is the analysis directory `/shared/myname/moringa`, and files end up here:

| Content | Location |
|---|---|
| Analysis results (00_Data … 07_Report) | `/shared/myname/moringa/results/` (`outdir: "results"`) |
| Script, logs and software versions of every step | `results/<step directory>/scripts/`, `logs/`, `logs/versions/` |
| Snakemake's own records, output of each SLURM job | `/shared/myname/moringa/.snakemake/` (job output in `.snakemake/slurm_logs/`) |
| Log of the main-process job | `/shared/myname/moringa/logs/` (set by `-o` in `run_annotation.sbatch`) |
| conda environments | `conda-prefix` in the profile |
| BRAKER4 images | `braker4.container_cache` |

The code directory `/shared/myname/annotation` is only read, never written; several species / analyses can share one
copy of the code. The directory from which `sbatch` is called only affects the relative path of the main-process log
(`-o logs/...`), not where the results go.

## 2. Directory layout

```
/shared/myname/
├── concordia/               workflow code (git clone --recursive)
├── db/                      databases (about 130 GB)
├── conda_envs/              conda environments (about 12 GB)
├── braker4_sif/             BRAKER4 images (about 3.4 GB)
└── moringa/                 one analysis directory per species (Snakemake's working directory)
    ├── config.yaml
    ├── rna_short.tsv  rna_long.tsv  homologs.tsv
    ├── slurm_profile/config.yaml     copied from concordia/profiles/slurm and edited
    ├── run_annotation.sbatch         script that submits the main process
    ├── logs/                         main-process logs
    └── results/                      created by the run
```

## 3. One-time setup (login node)

```bash
cd /shared/myname
# miniforge and the environment that runs the workflow
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/etc/profile.d/conda.sh
conda create -y -n concordia -c conda-forge -c bioconda \
    python=3.12 snakemake=8.18.2 snakemake-executor-plugin-slurm=2.6.0 pandas pyyaml pytest matplotlib graphviz

# workflow code (with the BRAKER4 submodule)
git clone --recursive <repository URL> concordia

# databases: one at a time, verified, can be re-run (resumes); InterProScan indexes are built on the first workflow run
bash concordia/scripts/download_databases.sh /shared/myname/db
```

Check that Apptainer is available on the compute nodes: `srun -p compute apptainer --version`. If the cluster uses
modules, it is `module load apptainer`; add that line to the script of step 7.

## 4. Configuration `moringa/config.yaml`

Copy `concordia/config/config.example.yaml` (Chinese comments: `config.example.zh.yaml`) and edit. Server-related settings:

```yaml
outdir: "results"                          # relative to the working directory: /shared/myname/moringa/results
species_prefix: "Mol"
genome: "/shared/myname/data/Moringa_oleifera.fasta"

rna:
  short: "rna_short.tsv"
  long:  "rna_long.tsv"
homologs: "homologs.tsv"
orthodb_proteins: "/shared/myname/db/orthodb/Viridiplantae.fa.gz"

repeat:
  dfam_dir: "/shared/myname/db/dfam"

braker4:
  container_cache: "/shared/myname/braker4_sif"
  bind: "/shared"                          # shared storage not mounted by Apptainer by default (must exist)
  executor: local                          # recommended: BRAKER4 runs inside the one job allocated to it (section 9)

eviann:
  swissprot: "/shared/myname/db/uniprot/uniprot_sprot.fasta"

qc:
  busco_lineages: [brassicales_odb12, embryophyta_odb12]
  busco_download_path: "/shared/myname/db/busco_downloads"   # shared with BRAKER4; missing datasets are downloaded (section 6)
  omark_db: "/shared/myname/db/omamer/LUCA.h5"

function:
  interproscan_dir: "/shared/myname/db/interproscan/interproscan-5.78-109.0"
  eggnog_data_dir: "/shared/myname/db/eggnog"
  go_obo: "/shared/myname/db/go/go-basic.obo"
  kegg_pathway_names: "/shared/myname/db/kegg/kegg_pathway_names.tsv"   # offline file when compute nodes have no internet

# Set according to the nodes: a single job cannot exceed one node (48 cores / 192 GB); runtime e.g. 12h / 3d
resources:
  default:      {threads: 4,  mem: 8G,    runtime: 12h}
  edta:         {threads: 48, mem: 128G,  runtime: 7d}
  braker4:      {threads: 48, mem: 160G,  runtime: 3d}
  eviann:       {threads: 32, mem: 64G,   runtime: 2d}
  gemoma:       {threads: 24, mem: 48G,   runtime: 2d}  # one job per homologous species, in parallel
  evm:          {threads: 32, mem: 32G,   runtime: 1d}
  repeatmasker: {threads: 32, mem: 32G,   runtime: 2d}
  interproscan: {threads: 16, mem: 32G,   runtime: 1d}  # function.chunks jobs in parallel
  eggnog:       {threads: 32, mem: 64G,   runtime: 1d}
  busco:        {threads: 16, mem: 32G,   runtime: 12h}
  omark:        {threads: 8,  mem: 16G,   runtime: 12h}  # LUCA.h5 needs about 12 GB
```

## 5. Profile `moringa/slurm_profile/config.yaml`

Copy `concordia/profiles/slurm/config.yaml` (Chinese comments: `config.zh.yaml`) and change these settings:

```yaml
executor: slurm
jobs: 50                                   # maximum jobs queued / running at the same time; follow cluster rules
default-resources:
  slurm_partition: "'compute'"             # partition name; keep the inner single quotes (see below)
  slurm_account: "'my_project'"            # when the cluster requires an accounting account; otherwise delete
  mem_mb: 8000                             # MB, integer only (units are not supported here)
  runtime: "10h"
software-deployment-method: [conda]
conda-prefix: /shared/myname/conda_envs    # absolute path on shared storage; set before the first build, never move it
conda-frontend: conda
latency-wait: 120
rerun-incomplete: true
restart-times: 1
```

**Quoting:** Snakemake evaluates every `default-resources` value as a Python expression. Text values such as
partition and account names therefore need single quotes inside the double quotes (`"'compute'"`): a bare name
that happens to be a Python built-in (`all`, `max`, `input` ...) would otherwise become a Python object, and every
job fails at start-up (e.g. `TypeError: all() takes no keyword arguments`). A dry-run (`-n`) does not catch this.

The resources requested by each step come from `resources` in `config.yaml`: `threads` → `--cpus-per-task`,
`mem` → `--mem`, `runtime` → `--time`. `-s` on the command line overrides the relative `snakefile` of the
profile, so the profile in the analysis directory needs no `snakefile` change.

## 6. Build environments and check (login node, needs internet)

```bash
source $HOME/miniforge3/etc/profile.d/conda.sh && conda activate concordia
cd /shared/myname/moringa
SM="snakemake -s /shared/myname/concordia/workflow/Snakefile -d . --configfile config.yaml --profile slurm_profile"

$SM --conda-create-envs-only      # builds the 7 environments (about 10-20 min); compute nodes then need no internet
conda clean -y --tarballs         # removes the downloaded archives (about 3 GB)

# OMArk downloads the NCBI taxonomy on its first run: fetch it in the omark environment (find it with $SM --list-conda-envs)
conda activate /shared/myname/conda_envs/<omark environment directory>
python -c "from ete3 import NCBITaxa; NCBITaxa()"
conda activate concordia

# BUSCO datasets of qc.busco_lineages missing from busco_download_path (or the files BRAKER4's compleasm needs):
# downloaded here, because the main process later runs as a job on a compute node that may have no internet
$SM busco_datasets

$SM -n                            # dry-run: checks the configuration and lists all jobs
```

## 7. Submit `moringa/run_annotation.sbatch`

The Snakemake main process only submits and tracks jobs and needs few resources, but it must run until the workflow
ends. Many clusters do not allow long processes on login nodes, so run it as a small job:

```bash
#!/bin/bash
#SBATCH -J annot_moringa
#SBATCH -p compute
#SBATCH -c 2
#SBATCH --mem 8G
#SBATCH -t 7-00:00:00            # not shorter than the whole workflow; if killed by the limit, resubmit to resume
#SBATCH -o logs/snakemake_%j.log
set -euo pipefail
# module load apptainer          # when the cluster uses modules
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate concordia
cd /shared/myname/moringa
snakemake -s /shared/myname/concordia/workflow/Snakefile -d . \
          --configfile config.yaml --profile slurm_profile
```

```bash
cd /shared/myname/moringa && mkdir -p logs
sbatch run_annotation.sbatch
squeue -u $USER                                  # the main process and the step jobs it submits
tail -f logs/snakemake_<jobid>.log               # progress
```

The step jobs inherit the environment of the main process (`PATH`, loaded modules ...), so conda and Apptainer only
need to be set up in this script.

## 8. During the run

- **Resuming**: after a failure or when the main process hit its time limit, `sbatch run_annotation.sbatch` again;
  completed steps are skipped.
- **Debugging a step**: read `results/<step directory>/logs/*.err` and `.snakemake/slurm_logs/`; to re-run it alone,
  get resources with `srun`, activate the environment named at the top of the script, and run
  `bash results/<step directory>/scripts/<script>.sh`.
- **Results**: final annotation `results/05_EVM/final/`, functional annotation `results/06_Function/merged/`
  (complete InterProScan result in `06_Function/interproscan/`, eggNOG in `06_Function/eggnog/`), report
  `results/07_Report/annotation_report_en.html` (Chinese: `annotation_report_zh.html`).

## 9. BRAKER4 `executor`

BRAKER4 runs nested: the workflow requests one job for it (here 48 cores / 160 GB) and starts another Snakemake in
it that runs BRAKER4's steps.

- **`executor: local` (recommended)**: BRAKER4's steps run inside that job, using its cores and memory; simplest, and
  independent of cluster policies.
- **`executor: slurm`**: BRAKER4's steps are submitted as SLURM jobs of their own. This needs `sbatch` on compute
  nodes (not allowed on many clusters), and the outer job still requests the full `resources.braker4` while only
  waiting for the sub-jobs, which wastes resources.
