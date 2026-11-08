# Concordia 示例：在 SLURM 集群上运行

[English](example_slurm.md) | 中文（两个版本内容一致，修改时请同步）

本示例假设：计算节点 48 核 / 192 GB 内存；登录节点可以联网，计算节点不一定；`/shared` 是所有节点都能访问的共享存储；
用户名以 `myname` 代替。路径、分区名和账户换成你们集群的即可。部署细节（数据库内容、常见问题）见 `deployment_guide.zh.md`。

## 1. 结果放在哪里

**结果写在 Snakemake 的工作目录里，而不是流程代码所在的目录。** 工作目录由 `-d` 指定；不加 `-d` 时就是运行
`snakemake` 时所在的目录。配置文件中的相对路径（包括 `outdir`、样本表、基因组）也都相对于工作目录解析；`~` 会展开为家目录，
但 `$USER` 等环境变量不会展开，请写实际路径。

本示例中，工作目录是分析目录 `/shared/myname/moringa`，各类文件的位置如下：

| 内容 | 位置 |
|---|---|
| 分析结果（00_Data … 07_Report） | `/shared/myname/moringa/results/`（`outdir: "results"`） |
| 每一步的脚本、日志、软件版本 | `results/<步骤目录>/scripts/`、`logs/`、`logs/versions/` |
| Snakemake 自己的记录、各 SLURM 作业的输出 | `/shared/myname/moringa/.snakemake/`（作业输出在 `.snakemake/slurm_logs/`） |
| 主进程作业的日志 | `/shared/myname/moringa/logs/`（由 `run_annotation.sbatch` 的 `-o` 指定） |
| conda 环境 | profile 中的 `conda-prefix` |
| BRAKER4 镜像 | `braker4.container_cache` |

流程代码目录 `/shared/myname/annotation` 只被读取，不写入任何结果，多个物种 / 多次分析可以共用同一份代码。
从哪个目录执行 `sbatch` 只影响主进程日志的相对路径（`-o logs/...`），不影响结果位置。

## 2. 目录规划

```
/shared/myname/
├── concordia/               流程代码（git clone --recursive）
├── db/                      数据库（约 130 GB）
├── conda_envs/              conda 环境（约 12 GB）
├── braker4_sif/             BRAKER4 镜像（约 3.4 GB）
└── moringa/                 一个物种一个分析目录（Snakemake 的工作目录）
    ├── config.yaml
    ├── rna_short.tsv  rna_long.tsv  homologs.tsv
    ├── slurm_profile/config.yaml     从 concordia/profiles/slurm 复制后修改
    ├── run_annotation.sbatch         提交主进程的脚本
    ├── logs/                         主进程日志
    └── results/                      运行后生成
```

## 3. 一次性部署（登录节点）

```bash
cd /shared/myname
# miniforge 与启动流程用的环境
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/etc/profile.d/conda.sh
conda create -y -n concordia -c conda-forge -c bioconda \
    python=3.12 snakemake=8.18.2 snakemake-executor-plugin-slurm=2.6.0 pandas pyyaml pytest matplotlib graphviz

# 流程代码（连同 BRAKER4 子模块）
git clone --recursive <仓库地址> concordia

# 数据库：逐个下载并校验，可重复运行（断点续传）；InterProScan 在流程第一次运行时自动建立索引
bash concordia/scripts/download_databases.sh /shared/myname/db
```

确认计算节点上有 Apptainer：`srun -p compute apptainer --version`。集群用 module 管理软件时为 `module load apptainer`，
并在第 7 步的脚本中加上这一句。

## 4. 配置文件 `moringa/config.yaml`

从 `concordia/config/config.example.zh.yaml`（中文注释；英文注释版为 `config.example.yaml`）复制后修改。与服务器有关的项：

```yaml
outdir: "results"                          # 相对于工作目录：/shared/myname/moringa/results
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
  bind: "/shared"                          # Apptainer 默认不挂载的共享存储（路径必须存在）
  executor: local                          # 推荐：BRAKER4 在分配给它的一个作业内运行（见第 9 节）

eviann:
  swissprot: "/shared/myname/db/uniprot/uniprot_sprot.fasta"

qc:
  busco_lineages: [brassicales_odb12, embryophyta_odb12]
  busco_download_path: "/shared/myname/db/busco_downloads"   # BRAKER4 也共用；缺少的数据集自动下载（第 6 节）
  omark_db: "/shared/myname/db/omamer/LUCA.h5"

function:
  interproscan_dir: "/shared/myname/db/interproscan/interproscan-5.78-109.0"
  eggnog_data_dir: "/shared/myname/db/eggnog"
  go_obo: "/shared/myname/db/go/go-basic.obo"
  kegg_pathway_names: "/shared/myname/db/kegg/kegg_pathway_names.tsv"   # 计算节点不能联网时用离线文件

# 按节点规格设置：单个作业不能超过一个节点（48 核 / 192 GB）；runtime 可写 12h / 3d 等
resources:
  default:      {threads: 4,  mem: 8G,    runtime: 12h}
  edta:         {threads: 48, mem: 128G,  runtime: 7d}
  braker4:      {threads: 48, mem: 160G,  runtime: 3d}
  eviann:       {threads: 32, mem: 64G,   runtime: 2d}
  gemoma:       {threads: 24, mem: 48G,   runtime: 2d}  # 每个同源物种一个作业，可并行
  evm:          {threads: 32, mem: 32G,   runtime: 1d}
  repeatmasker: {threads: 32, mem: 32G,   runtime: 2d}
  interproscan: {threads: 16, mem: 32G,   runtime: 1d}  # function.chunks 个作业并行
  eggnog:       {threads: 32, mem: 64G,   runtime: 1d}
  busco:        {threads: 16, mem: 32G,   runtime: 12h}
  omark:        {threads: 8,  mem: 16G,   runtime: 12h}  # LUCA.h5 约需 12 GB
```

## 5. Profile `moringa/slurm_profile/config.yaml`

复制 `concordia/profiles/slurm/config.zh.yaml`（中文注释；英文注释版为 `config.yaml`）为 `slurm_profile/config.yaml` 后修改以下几项：

```yaml
executor: slurm
jobs: 50                                   # 同时排队 / 运行的作业数上限，按集群规定调整
default-resources:
  slurm_partition: "'compute'"             # 分区名；保留里面的单引号（见下文）
  slurm_account: "'my_project'"            # 集群要求记账账户时填写；否则删掉这一行
  mem_mb: 8000
  runtime: 600
software-deployment-method: [conda]
conda-prefix: /shared/myname/conda_envs    # 共享存储上的绝对路径；第一次建环境前定好，之后不要移动
conda-frontend: conda
latency-wait: 120
rerun-incomplete: true
restart-times: 1
```

**引号：** Snakemake 把 `default-resources` 的每个值当作 Python 表达式求值。分区名、账户名这类文本值要在双引号里
再加单引号（`"'compute'"`）：否则恰好是 Python 内置名称的值（`all`、`max`、`input` 等）会变成 Python 对象，所有作业在
启动时失败（例如 `TypeError: all() takes no keyword arguments`）。dry-run（`-n`）检查不出这个问题。

各步骤作业申请的资源来自 `config.yaml` 的 `resources`：`threads` → `--cpus-per-task`，`mem` → `--mem`，`runtime` → `--time`。
命令行中的 `-s` 会覆盖 profile 里相对路径的 `snakefile`，所以分析目录中的 profile 不需要改 `snakefile`。

## 6. 建环境与检查（登录节点，需要联网）

```bash
source $HOME/miniforge3/etc/profile.d/conda.sh && conda activate concordia
cd /shared/myname/moringa
SM="snakemake -s /shared/myname/concordia/workflow/Snakefile -d . --configfile config.yaml --profile slurm_profile"

$SM --conda-create-envs-only      # 建好 7 个环境（约 10–20 分钟）；之后计算节点不需要联网
conda clean -y --tarballs         # 删除下载的压缩包（约 3 GB）

# OMArk 第一次运行要联网下载 NCBI 分类库：在 omark 环境中预先下载（环境目录可用 $SM --list-conda-envs 查到）
conda activate /shared/myname/conda_envs/<omark 环境目录>
python -c "from ete3 import NCBITaxa; NCBITaxa()"
conda activate concordia

# busco_download_path 中缺少的 qc.busco_lineages 数据集（及 BRAKER4 的 compleasm 需要的文件）在这里下载，
# 因为之后主进程作为作业运行在计算节点上，可能不能联网
$SM busco_datasets

$SM -n                            # dry-run：检查配置并列出全部任务
```

## 7. 提交运行 `moringa/run_annotation.sbatch`

Snakemake 主进程只负责提交和跟踪作业，占用资源很少，但要一直运行到流程结束。很多集群不允许在登录节点长时间运行进程，
所以用一个小作业运行它：

```bash
#!/bin/bash
#SBATCH -J annot_moringa
#SBATCH -p compute
#SBATCH -c 2
#SBATCH --mem 8G
#SBATCH -t 7-00:00:00            # 不能短于整个流程；超时被终止后重新提交即可续跑
#SBATCH -o logs/snakemake_%j.log
set -euo pipefail
# module load apptainer          # 集群用 module 管理时加上
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate concordia
cd /shared/myname/moringa
snakemake -s /shared/myname/concordia/workflow/Snakefile -d . \
          --configfile config.yaml --profile slurm_profile
```

```bash
cd /shared/myname/moringa && mkdir -p logs
sbatch run_annotation.sbatch
squeue -u $USER                                  # 主进程和它提交的各步骤作业
tail -f logs/snakemake_<jobid>.log               # 进度
```

各步骤作业继承主进程的环境（`PATH`、已加载的 module 等），所以 conda 和 Apptainer 只需在这个脚本中准备好。

## 8. 运行过程中

- **续跑**：中途失败或主进程超时后，重新 `sbatch run_annotation.sbatch`，已完成的步骤会跳过。
- **排查某一步**：看 `results/<步骤目录>/logs/*.err` 和 `.snakemake/slurm_logs/`；要单独重跑，用 `srun` 申请资源，
  激活脚本开头注明的环境后执行 `bash results/<步骤目录>/scripts/<脚本>.sh`。
- **结果**：最终注释 `results/05_EVM/final/`，功能注释 `results/06_Function/merged/`（完整 InterProScan 结果在 `06_Function/interproscan/`，eggNOG 在 `06_Function/eggnog/`），报告 `results/07_Report/annotation_report.html`。

## 9. BRAKER4 的 `executor`

BRAKER4 嵌套运行：流程为它申请一个作业（本例 48 核 / 160 GB），在其中再启动一个 Snakemake 运行 BRAKER4 的各步骤。

- **`executor: local`（推荐）**：BRAKER4 的步骤都在这个作业内运行，使用作业分到的核和内存；设置最简单，不受集群政策限制。
- **`executor: slurm`**：BRAKER4 的子步骤再各自提交为 SLURM 作业。要求计算节点上也能使用 `sbatch`（很多集群不允许）；
  而且外层作业仍按 `resources.braker4` 申请全部资源、却只在等待子作业，造成浪费。
