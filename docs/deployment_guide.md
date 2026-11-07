# Concordia 部署与数据库准备指南

[English](deployment_guide.en.md) | 中文（两个版本内容一致，修改时请同步）

> 适用版本：流程 v0.3.22 及之后（2026-11-06）。下载地址和软件版本均已于 2026-11-03 / 04 实际核实。

## 1. 资源规划

| | 全基因组（以 240 Mb 的植物基因组为例） | 小规模测试（单条染色体，约 13 Mb） |
|---|---|---|
| CPU | 32–48 核（SLURM 集群更好） | 12–14 核 |
| 内存 | 单节点 ≥ 128 GB（EDTA、BRAKER4）；OMArk 约 16 GB | 约 16 GB 可以运行（需相应调小 `resources`，OMArk 改用植物库 `Viridiplantae.h5`） |
| 数据库磁盘 | 约 130 GB（不含 FCS-GX）；FCS-GX 另需约 470 GB | 同左（数据库相同） |
| 结果磁盘 | 约 100–200 GB（BAM、EDTA 中间文件） | 约 10–20 GB |
| conda 环境 + 容器 | 约 25 GB | 同左 |
| 运行时间 | 2–5 天（集群并行） | 数小时到 1 天 |

磁盘建议：
- 数据库和结果放在大容量盘上；conda 环境放在本地 ext4/xfs 文件系统上（不要放在 NTFS 或网络盘）。
- 本流程的数据库、结果和容器缓存可以放在 NTFS 盘上，已验证可执行权限和符号链接正常。

## 2. 系统依赖

| 软件 | 用途 | 安装 |
|---|---|---|
| Linux x86-64 | — | — |
| conda / mamba（推荐 miniforge） | 所有生物信息软件 | https://github.com/conda-forge/miniforge |
| Apptainer（或 Singularity 3.x） | BRAKER4 的全部工具都在容器中运行 | 见下 |
| git、wget、curl | 获取代码和数据库 | 系统包管理器 |

### Apptainer

Ubuntu 22.04 / 24.04（需要 sudo）：

```bash
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
apptainer --version
```

注意：
- **Ubuntu 24.04** 默认限制非特权用户命名空间（`kernel.apparmor_restrict_unprivileged_userns=1`）。所以不能不用 root、直接从 conda 安装 Apptainer；请用上面的 PPA 安装，它会附带所需的 AppArmor 配置。
- **集群：** 通常已有 `module load apptainer` 或 `singularity`，直接使用即可。
- **验证：** `apptainer exec docker://alpine:3.20 echo ok` 应输出 `ok`。

## 3. 获取流程代码

```bash
# BRAKER4 是 git 子模块（external/BRAKER4），固定在已核实的 commit 6a5ee16，--recursive 会一并获取
git clone --recursive <本流程仓库> concordia && cd annotation
# 已经克隆但 external/BRAKER4 为空时，补拉子模块：
git submodule update --init
```

更新流程代码后，子模块的版本可能随之改变，执行 `git pull && git submodule update` 保持一致。

## 4. 运行环境

```bash
# 注意：conda/mamba 必须在 PATH 中（Snakemake 用它建立各步骤的环境），否则会报
#   "Command 'conda info --json' returned non-zero exit status 127"
export PATH=$HOME/miniforge3/bin:$PATH        # 或 conda init 后重新登录
# 运行流程本身用的环境（Snakemake 版本与 BRAKER4 要求一致）
mamba create -n concordia -c conda-forge -c bioconda \
    python=3.12 snakemake=8.18.2 snakemake-executor-plugin-slurm=2.6.0 pandas pyyaml pytest matplotlib
conda activate concordia

# 各步骤的软件环境由 Snakemake 按 workflow/envs/*.yaml 自动建立（版本已锁定）。
# 环境建在 profiles/slurm/config.yaml 的 conda-prefix 目录（默认 ~/annotation_conda_envs，约 15–20 GB）；
# 计算节点不能联网时，先在登录节点预先建好：
cd my_species && snakemake -s /path/to/concordia/workflow/Snakefile -d . --configfile config.yaml \
    --profile /path/to/concordia/profiles/slurm --conda-create-envs-only
```

`profiles/slurm/config.yaml` 中需要按集群修改的项：

| 项 | 说明 |
|---|---|
| `conda-prefix` | conda 环境存放目录。必须在共享文件系统上、登录节点和计算节点路径相同；环境按 yaml 内容和该目录真实路径的哈希命名，可被多个项目复用（移动该目录会导致环境全部重建；需要换位置时可用软链接）。默认 `~/annotation_conda_envs`（家目录在多数集群上是共享的）；家目录空间不足（约需 15–20 GB）时改到共享存储下 |
| `default-resources.slurm_partition` / `slurm_account` | 分区和记账账户 |
| `jobs` | 同时排队/运行的作业数上限 |

profile 已开启 `software-deployment-method: [conda]`，命令行不必再加 `--sdm conda`。

配置模板和每个 profile 都有英文注释版（`config/config.example.en.yaml`、`profiles/*/config.en.yaml`），参数与中文版完全相同（`tests/test_config_translations.py` 检查两者一致）。Snakemake 只读取 profile 目录下的 `config.yaml`，想用英文版时把 `config.en.yaml` 复制为 `config.yaml`。

环境清单和版本见 `workflow/envs/README.md`。其中：
- 默认建立 7 个环境（core、repeat、eviann、agat、qc、omark、function）；repeat、eviann、agat、qc 与其它软件存在依赖或文件冲突，必须独立（原因见 `workflow/envs/README.md`）。
- InterProScan 需要 Java 11+，由 `function` 环境中的 openjdk 11 提供，不依赖系统 Java。

## 5. 数据库

### 5.1 一键下载

```bash
bash scripts/download_databases.sh /path/to/databases --limit-rate 30m   # 可重复运行，断点续传
```

脚本**逐个**处理：一个文件下载、校验、解压、`sync` 完成后才开始下一个；每步开始前检查内核日志，如果 USB 存储报告命令超时就暂停等待。之所以不并行，是因为在 USB 移动硬盘上同时下载和解压多个大文件，可能使部分 USB 转接芯片（UAS 驱动下）反复超时、硬盘停止响应（见第 8 节）。完成后的目录结构如下：

| 目录 | 内容 | 下载大小 | 解压后 | 配置项 |
|---|---|---|---|---|
| `orthodb/Viridiplantae.fa.gz` | OrthoDB v12 植物蛋白（BRAKER 分区；其它类群换成相应分区，如 Metazoa、Fungi） | 3.1 GB | 不解压 | `orthodb_proteins` |
| `orthodb/odb12v2_levels.tab.gz`、`odb12v2_level2species.tab.gz` | OrthoDB 分类阶元表：按阶元（如 Brassicales）抽取蛋白子集，减少 BRAKER4 的磁盘占用 | 0.3 MB | — | `orthodb_subset.levels_tab`、`orthodb_subset.level2species_tab` |
| `omamer/LUCA.h5` | OMAmer 数据库（OMArk 用），运行约需 12 GB 内存 | 9.9 GB | — | `qc.omark_db` |
| `omamer/Viridiplantae.h5` | 植物 OMAmer 数据库：内存不足 16 GB 时代替 LUCA.h5（约需 2 GB 内存；只能识别植物来源的污染） | 1.5 GB | — | `qc.omark_db` |
| `busco_downloads/lineages/<lineage>` | BUSCO 数据集（离线模式）；BRAKER4 的 compleasm 也使用这个目录 | 每个约 0.2 GB | — | `qc.busco_download_path`、`qc.busco_lineages` |
| `busco_downloads/lineages/` 中的 `eukaryota_odb12`、`placement_files/`、`file_versions.tsv`、`*.done` | BRAKER4 的 compleasm 额外需要的文件（`compleasm` 步骤下载；目录结构与 compleasm 自己下载的一致） | 约 0.2 GB | — | 同上 |
| `eggnog/` | eggNOG-mapper 5.0.2：`eggnog.db`、`eggnog_proteins.dmnd`、taxa | 12 GB | 约 58 GB | `function.eggnog_data_dir` |
| `interproscan/interproscan-5.78-109.0/` | InterProScan 完整发行版（含全部成员数据库） | 7.0 GB | 约 30–60 GB | `function.interproscan_dir` |
| `dfam/dfam40.0.h5` | Dfam 4.0 FamDB 根分区：RepeatMasker 4.2 即使用自定义库（EDTA）也需要；流程会自动把它链接进 conda 环境 | 60 MB | 450 MB | `repeat.dfam_dir` |
| `go/go-basic.obo` | GO 本体（报告中的 GO 名称和 BP/MF/CC 分类） | 32 MB | — | `function.go_obo` |
| `uniprot/uniprot_sprot.fasta` | UniProt-SwissProt：EviAnn 只用 RNA 数据注释时用于确定 CDS（不填则 EviAnn 改用同源物种蛋白） | 90 MB | 约 280 MB | `eviann.swissprot` |
| `kegg/kegg_pathway_names.tsv` | KEGG 通路名称（KEGG REST，仅限学术用途） | 很小 | — | `function.kegg_pathway_names` |

### 5.2 各数据库说明

- **BRAKER4 与 BUSCO 数据集：** BRAKER4 在基因预测中运行 compleasm（BUSCO 基因的 CDS 提示给 AUGUSTUS；TSEBRA 合并后补回丢失的 BUSCO 基因），这与流程的 BUSCO 评估无关、不能关闭（`skip_busco` 只关闭 BRAKER4 自己的 BUSCO 评估）。流程把 `qc.busco_download_path` 写入 BRAKER4 的 `config.ini`，两者共用同一份数据集，并自动挂载进 BRAKER4 的容器；compleasm 每次运行会刷新该目录下的 `file_versions.tsv.hash`，所以目录需要可写。未设置 `qc.busco_download_path` 时，BRAKER4 下载到结果目录 `03_Ab_initio/BRAKER4/busco_downloads`（不再写入 `external/BRAKER4`）。
- **BUSCO lineage 的选择：** 先在 https://busco-data.ezlab.org/v5/data/file_versions.tsv 里确认名称存在。
  - 例如辣木属于十字花目，应选 `brassicales_odb12`，再加更广的 `embryophyta_odb12`。
  - **odb12 中没有 `eudicots_odb12`**，只有 `eudicots_odb10`。
  - 下载脚本中的 lineage 列表可按物种修改。
- **InterProScan：**
  - 解压后需要运行一次 `python3 setup.py -f interproscan.properties` 来建立 HMM 索引。流程中的 `interproscan_setup` rule 会在第一次运行时自动完成，并在安装目录写入 `.pipeline_setup_done` 标记。
  - 需要 Java 11+。
  - 完整发行版包含全部成员数据库；需要时可以用 `function.interproscan_appl` 只运行部分数据库，例如 `Pfam,PANTHER`。
- **eggNOG-mapper：** 也可以用 `download_eggnog_data.py --data_dir <dir> -y` 下载；脚本直接从 `http://eggnog5.embl.de/download/emapperdb-5.0.2/` 下载相同的文件。
- **OrthoDB：** BRAKER4 推荐用目标物种所在大类群的分区（植物用 Viridiplantae），流程会把它与同源物种蛋白合并、去重。
- **KEGG 通路名：** 计算节点不能联网时，填写已下载文件的路径；设为 `download` 时会在运行时下载。
- **NCBI FCS-GX（可选，默认关闭）：**
  - 数据库约 470 GB，运行时建议放在 `/dev/shm`，需要约 500 GB 内存。
  - 下载命令（依据 NCBI FCS 文档，安装 bioconda `ncbi-fcs-gx` 后）：
    ```bash
    sync_files.py get --mft https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest --dir /path/to/gxdb
    ```
  - 在配置中填写 `assembly_qc.fcs_gx.db`（例如 `/path/to/gxdb/all`）和 `tax_id`。

### 5.3 在线服务（不推荐用于整个基因组）

| 工具 | 在线方式 | 限制 |
|---|---|---|
| InterProScan | EBI Job Dispatcher REST 接口 | 必须提供邮箱；同时提交不超过 30 个任务；逐条蛋白提交，3 万个蛋白需要数天 |
| eggNOG-mapper | 官方网页服务器 | 只能网页提交，结果通过邮件返回，无法嵌入流程 |

## 6. 配置

```bash
mkdir my_species && cd my_species           # 分析目录：配置、样本表和结果都在这里，代码目录只被读取
cp /path/to/concordia/config/config.example.yaml config.yaml
```

必填项：

| 项 | 说明 |
|---|---|
| `genome`、`species_prefix` | 基因组和基因 ID 前缀 |
| `rna.short` / `rna.long` | 样本表（长读须注明 `platform`：`pacbio_isoseq` / `pacbio_hifi` / `ont_cdna` / `ont_drna`） |
| `homologs` | 同源物种的基因组、GFF3、蛋白 |
| `orthodb_proteins` | OrthoDB 分区蛋白 |
| `qc.busco_lineages`、`qc.busco_download_path`、`qc.omark_db` | 质量评估 |
| `function.interproscan_dir`、`function.eggnog_data_dir` | 功能注释 |
| `braker4.executor` | 在集群上设为 `slurm`，BRAKER4 的子任务也会提交到 SLURM |
| `braker4.bind` | 容器需要访问的路径（数据库、结果目录所在的挂载点） |

示例：仓库中的 `tests/config.test.yaml`（模拟数据，用于 dry-run）；分析目录的组织方式见 README"快速开始"。

## 7. 运行

```bash
conda activate concordia
cd my_species; A=/path/to/concordia
SM="snakemake -s $A/workflow/Snakefile -d . --configfile config.yaml"
# 1) 检查
$SM --profile $A/profiles/slurm -n
# 2) 运行（SLURM；完整示例见 example_slurm.md）
$SM --profile $A/profiles/slurm
# 单机
$SM --profile $A/profiles/local       # 环境同样建在 ~/annotation_conda_envs，可与 SLURM 共用
# 只跑到某一步
$SM --profile $A/profiles/local --until evm
$SM --profile $A/profiles/local data_qc
```

- **断点续跑：** 重新运行同一命令即可。
- **排查某一步：** 看 `<步骤目录>/logs/` 中的 `.err`，再激活脚本头部写明的环境，`bash <步骤目录>/scripts/<编号>_<rule>.sh` 单独重跑。
- **完整示例：** 在 SLURM 集群上从部署到提交的完整步骤见 `example_slurm.md`。
- **结果：** `<outdir>/05_EVM/final/` 是最终注释，`<outdir>/07_Report/annotation_report.html` 是统计报告。

## 8. 常见问题

| 现象 | 原因与处理 |
|---|---|
| BRAKER4 报 `apptainer: command not found` 或无法创建命名空间 | 没有安装 Apptainer，或者在 Ubuntu 24.04 上从 conda 安装（受 AppArmor 限制）。按第 2 节用 PPA 安装 |
| GeMoMa 只打印工具列表就退出 | 在 bioconda 包装脚本后面又写了 `CLI`；`gemoma.cmd` 应为 `GeMoMa -Xmx{mem_gb}g` |
| BUSCO 报找不到 lineage | lineage 名称不存在（例如 `eudicots_odb12`），或者 `busco_download_path` 没有指向包含 `lineages/` 的目录 |
| OMArk 被系统 OOM 终止 | OMAmer 会把整个 LUCA.h5 读进内存；`resources.omark.mem_mb` 至少 13 GB，并保证它不和其它大任务同时运行 |
| InterProScan 报 Java 版本过低 | 系统 Java 是 8；流程使用 `function` 环境里的 openjdk 11，请确认 profile 中启用了 `software-deployment-method: [conda]`（两个 profile 默认已启用），没有改用 PATH 中的软件 |
| 建环境时报 `create: The following argument was not expected: --no-default-packages` | Snakemake 8.18.2 与 mamba 2.x 不兼容；profile 中已设置 `conda-frontend: conda`（conda 默认使用 libmamba 求解器） |
| BRAKER4 运行中把系统盘写满 | BRAKER4 会把蛋白库复制约 3 份（运行目录、GeneMark-ETP 的 arx/ 和 data/）并生成 ProtHint 临时文件：OrthoDB 植物全库（5.4 GB）约需 20 GB。磁盘有限时用 `orthodb_subset.levels` 只取近缘阶元（例如 Brassicales + Malvales 约 206 万条） |
| OMArk（omamer search）被 OOM 杀掉（`已杀死` / `Out of memory: Killed process ... (omamer)`） | OMAmer 把整个数据库读入内存，LUCA.h5 约需 12 GB。内存小的机器改用 `Viridiplantae.h5`，或把 `resources.omark.mem_mb` 设到 14000 以上，使它不与其它任务同时运行 |
| AGAT 报 `Cannot open agat_tmp_<文件名>/result_*.stor` 或 `Failed to delete agat_tmp_<文件名>` | AGAT 在当前目录下按输入文件名建临时文件夹，两个任务同时处理同一个文件时互相删除对方的临时文件。流程中每个 AGAT 任务都在自己的临时目录（`<脚本>.agat`）中运行 |
| InterProScan 报 `'sequence' is not an amino acid sequence [...]` | 蛋白中有内部终止密码子（gffread 写成 `.`）。流程把功能注释输入中的内部终止写为 X，并在最终 GFF3 中把这类转录本标记为 `internal_stop=yes`（`05_EVM/final/internal_stops.tsv`） |
| OMArk 第一次运行时卡在 `Downloading taxdump.tar.gz from NCBI FTP site`，或在不能联网的节点上失败 | ete3 第一次使用时下载 NCBI 分类库（约 80 MB 下载、780 MB 的 `~/.etetoolkit/taxa.sqlite`）。先在能联网的节点上、omark 环境中运行一次 `python -c "from ete3 import NCBITaxa; NCBITaxa()"` |
| 报错 `at least one of rna.short / rna.long must list samples` | 没有找到任何 RNA 样本：检查 `rna.short` / `rna.long` 的样本表路径和内容。确实没有 RNA-seq 时设 `rna.allow_none: true`（无 RNA 模式，见配置模板中 `rna.allow_none` 的说明） |
| AGAT 报 `progress bar already finished at .../AGAT/OmniscientI.pm line 390` | AGAT 1.7.0 的进度条计数超过预估行数后直接退出（GeMoMa 输出上可复现）。流程生成 `00_Data/agat_config.yaml`（`progress_bar: false`），所有 AGAT 调用都通过 `--config` 使用它 |
| EDTA 报 `Error: The RMblast engine is not installed in RepeatMasker!` | 误导性的提示：真正原因是 RepeatMasker 4.2 找不到 Dfam FamDB（`FamDB data directory not found`），bioconda 包不带这份数据。设置 `repeat.dfam_dir`，`repeat_famdb_setup` 步骤会把它链接到环境里 famdb 的 `Libraries/famdb/` |
| EDTA 报 `Raw TIR results not found`，`.err` 中有 TIR-Learner 的 `KeyError: 0`（`check_TIR_TSD.py`） | TIR-Learner 3.0.7 不兼容 pandas 3（`Series[0]` 不再按位置取值）；`repeat` 环境已锁定 `pandas=2.3.3` |
| `calcDivergenceFromAlign.pl` 报 `Can't locate RepeatMaskerConfig.pm` | bioconda 把 RepeatMasker 的工具脚本复制到 `bin/`，脚本按自身上级目录找模块；rule 中已将 `PERL5LIB` 指向 RepeatMasker 的实际安装目录 |
| `conda info --json` 返回 127 | conda 不在 PATH 中；执行 `export PATH=$HOME/miniforge3/bin:$PATH` 或 `conda init` |
| USB 移动硬盘在大量读写时卡死，内核日志反复出现 `uas_eh_abort_handler`，进程停在 `D` 状态 | 部分 USB 转接芯片在 UAS 驱动下不稳定。处理：① 改用 usb-storage 驱动：用 `lsusb` 查到硬盘盒的 `<VID>:<PID>`，`echo "options usb-storage quirks=<VID>:<PID>:u" \| sudo tee /etc/modprobe.d/disable-uas.conf && sudo update-initramfs -u`，重启或重新插拔后生效；② 直接接电脑的 USB 口，不经过扩展坞；③ 移动硬盘只存放数据库，数据和结果放在内置硬盘；④ 不要在移动硬盘上并行解压多个大文件 |
| GeneMark 许可 | GeneMark 对商业用途有许可限制，商业项目请先确认许可 |
