# Concordia：真核基因组结构与功能注释流程（Snakemake）

[English](README.md) | 中文

Concordia（拉丁语"一致"）把 RNA、从头预测和同源投射三类证据整合成共识注释。输入目标基因组、RNA-seq（短读和/或 PacBio / ONT 长读）以及若干同源物种的基因组和注释，依次完成以下步骤，最后生成一份图文统计报告：
- 组装质量评估
- 重复序列注释与软屏蔽
- RNA 比对
- 三路基因预测（EviAnn / BRAKER4 / GeMoMa）
- EVM 整合
- 与 EviAnn 比较修正
- 过滤
- 功能注释

- 部署指南：`docs/deployment_guide.zh.md`（英文版 `docs/deployment_guide.md`）
- SLURM 集群运行示例：`docs/example_slurm.zh.md`（英文版 `docs/example_slurm.md`）
- 中英文文档（README、部署指南、SLURM 示例、配置模板、profile）内容一致，修改时请同步修改两个版本。

## 流程概览

| 目录 | 内容 | 主要工具 |
|---|---|---|
| `00_Data` | 输入汇集（可按 OrthoDB 阶元抽取 BRAKER4 用的蛋白子集，`orthodb_subset`；未改动的输入用软链接，不复制：基因组、单个同源蛋白/CDS 文件、只用 OrthoDB 时的蛋白库）；序列名超过 13 个字符时统一重命名；同源注释只保留蛋白编码基因；组装质量评估（连续性、BUSCO、端粒、LAI（EDTA 完成后计算），可选 Merqury / FCS-GX）；RNA 原始数据质控 | seqkit, BUSCO, tidk, fastp, NanoPlot |
| `01_Repeat` | TE 库构建与软屏蔽；TE 组成、Kimura 景观、LTR 插入时间统计 | EDTA, RepeatMasker |
| `02_RNA` | 短读和长读比对、合并 BAM、StringTie 组装、长读内含子链；EviAnn 注释（只用 RNA 数据；CDS 判定用 SwissProt，未配置时退回同源蛋白） | HISAT2, minimap2, samtools, StringTie2, EviAnn |
| `03_Ab_initio` | BRAKER4（嵌套 Snakemake）：短读和 / 或 PacBio 长读加蛋白证据，模式按数据自动选择（dual / ETP / IsoSeq / EP，见下文）；ONT 长读不进入 BRAKER4 | BRAKER4 |
| `04_Homology` | 同源投射与蛋白比对：每个同源物种单独运行 GeMoMa 和 miniprot（各自做 BUSCO/OMArk），再合并（GeMoMa GAF；miniprot 拼接）进入 EVM | GeMoMa, miniprot |
| `05_EVM` | EVM 整合 → 修正（R1–R6）→ TE / 长度过滤 → 按染色体命名（如 `Mol01G000010`）→ 还原序列名 → 导出；各来源的 BUSCO / OMArk / AGAT 评估 | EVidenceModeler 2.1（Rust 版可选）, AGAT, TEsorter |
| `06_Function` | 结构域、GO、KEGG、COG（完整 InterProScan 结果 `interproscan/<前缀>.interproscan.{tsv,gff3}`、eggNOG `eggnog/`、合并表 `merged/`） | InterProScan, eggNOG-mapper |
| `07_Report` | `annotation_report_zh.html`（中文）和 `annotation_report_en.html`（英文）：每一步都有统计表和图，内容相同 | Python |

**流程图：**
- `docs/workflow_main.png` / `.svg`：主流程（有 RNA-seq），按 00–07 目录分组，不含质控和统计步骤；箭头表示"读取其输出"，所有直接依赖都画出。
- `docs/workflow_main_norna.png` / `.svg`：无 RNA 模式（`rna.allow_none: true`）的主流程。
- `docs/workflow_rulegraph.png` / `.svg`：全部步骤（示例配置下 103 个，质控和统计步骤为浅灰色）。
- 修改规则后重新生成：`bash scripts/make_workflow_diagrams.sh`（需要 snakemake 和 graphviz；用 `tests/config.test.yaml` 示例配置，不需要数据库）。

**每一步的命令都会先写成脚本：** `<目录>/scripts/<编号>_<rule>[.<通配符>].sh`，然后再执行；日志写入 `<目录>/logs/`，软件版本写入 `<目录>/logs/versions/`，运行时间、CPU 时间和峰值内存（Snakemake benchmark）写入 `<目录>/logs/benchmarks/`，报告最后一节将其与申请的资源对比并给出调整建议。排查问题时，激活脚本头部注明的环境，直接 `bash` 该脚本即可单独重跑这一步。

## 快速开始

结果写在 Snakemake 的工作目录（`-d`）里，不写进流程代码目录。建议流程代码只放一份，每个物种 / 每次分析各用一个目录。

```bash
# 1. 获取流程代码（BRAKER4 是 git 子模块，固定在 commit 6a5ee16，--recursive 一并获取）
git clone --recursive <本流程仓库> concordia      # 已克隆但缺少 BRAKER4 时：git submodule update --init
A=$PWD/concordia                                  # 流程代码目录，下面用到
#   （可选）git clone https://github.com/hyphaltip/EVidenceModeler_rust $A/external/EVidenceModeler_rust

# 2. 新建分析目录，放配置和样本表
mkdir my_species && cd my_species
cp $A/config/config.example.zh.yaml config.yaml    # 中文注释版；英文注释版为 config.example.yaml（参数相同）
cp $A/config/rna_short.example.tsv rna_short.tsv
cp $A/config/rna_long.example.tsv  rna_long.tsv
cp $A/config/homologs.example.tsv  homologs.tsv
#    编辑 config.yaml：genome、样本表（模板默认即 rna_short.tsv、rna_long.tsv、homologs.tsv）、
#    qc.busco_lineages、qc.omark_db、function.interproscan_dir、function.eggnog_data_dir；相对路径都相对于本目录

# 3. 运行（-d . 让结果写在本目录的 results/ 下）
SM="snakemake -s $A/workflow/Snakefile -d . --configfile config.yaml"
$SM --profile $A/profiles/local -n                 # dry-run
$SM --profile $A/profiles/local                    # 单机；按 workflow/envs 建 conda 环境（~/annotation_conda_envs）
$SM --profile $A/profiles/slurm                    # SLURM；完整步骤见 docs/example_slurm.zh.md

# 只跑某一部分
$SM --profile $A/profiles/local data_qc            # 00_Data 评估
$SM --profile $A/profiles/local --until evm        # 到 EVM 为止
```

软件清单与版本见 `workflow/envs/README.md`，各环境的版本均已锁定。

### 分析目录的内容

运行后，分析目录中的配置里的相对路径、`results/` 和 `.snakemake/`（Snakemake 自己的记录）都在分析目录内，流程代码目录只被读取：

```
my_species/
├── config.yaml
├── rna_short.tsv  rna_long.tsv  homologs.tsv
├── .snakemake/          Snakemake 的运行记录（断点续跑依赖它）
└── results/             00_Data ... 07_Report
```

## 根据数据自动做出的选择

| 数据 | BRAKER4 | EviAnn 输入（每个样本一行） | EVM 转录本证据 | GeMoMa RNA 证据 | 修正步骤的 RNA 支持 |
|---|---|---|---|---|---|
| 短读 + PacBio (+ONT) | dual | 短读：bam；长读：bam_isoseq | 长读 StringTie | 短读 BAM | 长读内含子链 |
| 短读 + ONT | ETP（不用 ONT） | 短读：bam；ONT：bam_isoseq（暂定） | 长读 StringTie | 短读 BAM | 长读内含子链 |
| 只有短读 | ETP | bam | 短读 StringTie | 短读 BAM | 短读：剪接 reads（B）+ StringTie 内含子链（A）+ 单外显子（C） |
| 只有 PacBio | IsoSeq | bam_isoseq | 长读 StringTie | 无（可选长读） | 长读内含子链 |
| 只有 ONT | EP | bam_isoseq（暂定） | 长读 StringTie | 无（可选长读） | 长读内含子链 |
| 没有 RNA（需 `rna.allow_none: true`） | EP | 不运行 | 无 | 无 | 无：修正只做 R5 证据检查；过滤以蛋白证据代替 RNA 证据；无 UTR / 可变剪接 / lncRNA |

EviAnn 拿到的是每个样本各自的比对结果，每个样本一行（不是合并后的 BAM）：它对每个样本分别组装，lncRNA 只有在
至少 2 个样本中检测到（TPM ≥ 0.5）才会保留，编码转录本要在至少 1/10 的样本中检测到（TPM ≥ 100 时除外）；全部 RNA
只有 1 个样本时不会输出 lncRNA。样本多为少数组织的重复时，填写样本表中可选的 `tissue` 列（或 `sample_group`）：
同一组织的样本合并后作为一个样本交给 EviAnn，同一组织的短读和长读一起给出（`bam_mix`）。只影响 EviAnn 的输入。

这些选择的逻辑在 `workflow/lib/common.py`，测试在 `tests/test_common.py`。

## 关键规则

**EVM 权重调优（`05_EVM/tune`，配置 `evm.tune`，默认关闭）**
- 默认关闭：测试中最优权重的 F 只比配置的权重高 0.001–0.003，却要多花约 11 CPU 小时；`evm.tune.enabled: true` 打开。
- 在一条序列（默认最长者）上，对各证据来源的候选权重做网格搜索，默认 16 组加配置中的基线。
- 用原始证据给每组结果打分：P 为多外显子模型的内含子链被长读 / 短读（与修正步骤同一标准）/ miniprot 任一完全支持的比例，R 为证据位点被覆盖的比例，F = 2PR/(P+R)。
- 默认用 F 最高的权重跑全基因组 EVM。

**修正（`05_EVM/refine`，`workflow/scripts/refine.py`）**
- 按 CDS 内含子链比较 EVM 和 EviAnn 模型。
- 长读支持以"全长 reads"计：Iso-Seq HQ 转录本按 `rna_long.tsv` 的 `counts` 列（isoseq cluster 的 `cluster_report.csv` 或 isoseq collapse 的 `flnc_count.txt`）或表头的 `full_length_coverage` 加权，原始 reads（CCS/FLNC、ONT）按 1 计（`refine.lr_weight_by_coverage`）。长读可以是 FASTA 或 FASTQ。HQ 转录本既无计数文件又无表头标签时支持度会被低估，报告中会提示；此时可将 `refine.min_lr_support` 设为 1。
- 规则：R1 引入 UTR 和异构体；R2 补入缺失基因；R3 冲突时用有长读支持的模型替换；R4 处理拆分和融合；R5 检查 EVM 独有基因的证据；R6 lncRNA 单独输出。
- 每个位点的处理写入 `refine_decisions.tsv`。需要人工检查的位点会在报告中单独列出。

**命名（`05_EVM/final`，`workflow/scripts/rename_genes.py`，配置 `naming`）**
- 编码基因 `<前缀><染色体号>G<编号>`，如 `Mol01G000010`；lncRNA 用 `L`，如 `Mol01L000010`。
- 每条染色体内按位置编号，步长 10。
- 染色体号从原始序列名中取（Chr1 / chr01 / Chromosome_3 / 1）；匹配不到的 scaffold 统一用 `00`；也可以手动指定，如 `ChrM: Mt`。
- 转录本为 `基因.1`、`.2`（`.1` 为 CDS 最长者）。
- 主结果 `<前缀>.gff3`（及 cds/pep/cdna）每个基因只保留一个转录本（CDS 最长，等长时取转录本最长）；全部异构体另存为 `<前缀>.all_isoforms.*`。最终结果的统计和评估都基于主结果。
- 补出 5′/3′ UTR 特征；新旧 ID 对照表为 `id_map.tsv`。

**过滤（`05_EVM/final`）**
- 基因 CDS 与 TE 重叠超过 `max_te_overlap`（默认 0.5），**并且**没有 RNA 证据，才会被剔除。
- RNA 证据满足任一条即可：StringTie、来自 EviAnn、BRAKER4 gene_support。
- CDS 短于 `min_cds_aa` 的模型单独剔除；但有 RNA 证据（无 RNA-seq 时为蛋白证据）的保留，并标记 `short_cds=yes`（`refine.keep_short_supported`，默认开启）。
- TEsorter 结果只记录，不作为删除条件。
- 每个基因的判断写入 `filter_decisions.tsv`。

## 测试

```bash
python -m pytest tests                                                     # 单元测试
snakemake -s workflow/Snakefile --configfile tests/config.test.yaml -n     # 假数据 dry-run
```

## 许可证

Concordia 的代码以 [MIT 许可证](LICENSE) 发布（Copyright (c) 2026 Jiyang Chang）。

### 第三方软件与数据库

Concordia 只通过命令行调用下列软件，不包含、也不分发它们的代码；它们由 conda 或各自的容器 / 安装包安装，
各自遵守自己的许可证。使用本流程即同时使用这些软件，发表结果时请一并引用。商业用途请特别注意标 **⚠** 的项。

| 软件 | 许可证（据 bioconda 包信息） |
|---|---|
| BRAKER4（子模块，`external/BRAKER4`） | MIT |
| **⚠ GeneMark**（GeneMark-ETP / ES，在 BRAKER4 容器内） | 学术免费；**商业用途需向作者另行取得许可** |
| AUGUSTUS、ProtHint、TSEBRA、compleasm 等（BRAKER4 容器内） | 各自的开源许可证，见 BRAKER4 文档 |
| EDTA 2.3.0、LTR_retriever、TIR-Learner、TEsorter | GPL-3.0 |
| RepeatMasker 4.2.4、RepeatModeler 2.0.9 | Open Software License 2.1 |
| EviAnn 2.0.7 | GPL-3.0 |
| HISAT2 2.2.3 | GPL-3.0 |
| minimap2、samtools、StringTie、miniprot、gffread、seqkit、fastp、NanoPlot、tidk、BUSCO | MIT |
| gffcompare | Artistic License 2.0 |
| GeMoMa 1.9 | GPL-3.0 |
| EVidenceModeler 2.1.0 | BSD-3-Clause |
| AGAT 1.7.0 | GPL-3.0 |
| OMArk、OMAmer | LGPL-3.0 |
| **⚠ eggNOG-mapper 2.1.12** | AGPL-3.0（以网络服务形式提供时有额外的开源义务） |
| diamond、MultiQC | GPL-3.0 |
| InterProScan | Apache-2.0；部分成员数据库 / 工具（如 SignalP、TMHMM、Phobius）需另行取得许可，默认发行版不含 |
| NCBI FCS-GX | 公共领域（NCBI） |
| Snakemake | MIT |

数据库：UniProt-SwissProt 为 CC BY 4.0，Dfam 为 CC0；**⚠ KEGG REST 接口仅限学术用途**（`function.kegg_pathway_names`）；
OrthoDB、OMA、eggNOG、BUSCO 数据集、InterPro、GO 等请遵守各自网站上的使用条款。
上表根据 2026-11 安装的版本整理，仅供参考，以各软件 / 数据库的官方说明为准。
