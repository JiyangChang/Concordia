# 真核基因组结构与功能注释流程（Snakemake）

[English](README.en.md) | 中文

输入目标基因组、RNA-seq（短读和/或 PacBio / ONT 长读）以及若干同源物种的基因组和注释，依次完成以下步骤，最后生成一份图文统计报告：
- 组装质量评估
- 重复序列注释与软屏蔽
- RNA 比对
- 三路基因预测（EviAnn / BRAKER4 / GeMoMa）
- EVM 整合
- 与 EviAnn 比较修正
- 过滤
- 功能注释

- 部署指南：`docs/deployment_guide.md`（英文版 `docs/deployment_guide.en.md`）
- 中英文文档（README、部署指南、配置模板、profile）内容一致，修改时请同步修改两个版本。

## 流程概览

| 目录 | 内容 | 主要工具 |
|---|---|---|
| `00_Data` | 输入汇集（可按 OrthoDB 阶元抽取 BRAKER4 用的蛋白子集，`orthodb_subset`；未改动的输入用软链接，不复制：基因组、单个同源蛋白/CDS 文件、只用 OrthoDB 时的蛋白库）；序列名超过 13 个字符时统一重命名；同源注释只保留蛋白编码基因；组装质量评估（连续性、BUSCO、端粒、LAI（EDTA 完成后计算），可选 Merqury / FCS-GX）；RNA 原始数据质控 | seqkit, BUSCO, tidk, fastp, NanoPlot |
| `01_Repeat` | TE 库构建与软屏蔽；TE 组成、Kimura 景观、LTR 插入时间统计 | EDTA, RepeatMasker |
| `02_RNA` | 短读和长读比对、合并 BAM、StringTie 组装、长读内含子链；EviAnn 注释（只用 RNA 数据；CDS 判定用 SwissProt，未配置时退回同源蛋白） | HISAT2, minimap2, samtools, StringTie2, EviAnn |
| `03_Ab_initio` | BRAKER4（嵌套 Snakemake；只使用 PacBio 长读） | BRAKER4 |
| `04_Homology` | 同源投射与蛋白比对：每个同源物种单独运行 GeMoMa 和 miniprot（各自做 BUSCO/OMArk），再合并（GeMoMa GAF；miniprot 拼接）进入 EVM | GeMoMa, miniprot |
| `05_EVM` | EVM 整合 → 修正（R1–R6）→ TE / 长度过滤 → 按染色体命名（如 `Mol01G000010`）→ 还原序列名 → 导出；各来源的 BUSCO / OMArk / AGAT 评估 | EVidenceModeler 2.1（Rust 版可选）, AGAT, TEsorter |
| `06_Function` | 结构域、GO、KEGG、COG | InterProScan, eggNOG-mapper |
| `07_Report` | `annotation_report.html`：每一步都有统计表和图 | Python |

**流程图：**
- `docs/workflow_main.png` / `.svg`：主流程，按 00–07 目录分组，不含质控和统计步骤。
- `docs/workflow_rulegraph.png` / `.svg`：全部 90 个步骤。
- 重新生成：`snakemake -s workflow/Snakefile --configfile tests/config.test.yaml --rulegraph | dot -Tsvg > rulegraph.svg`（需要 graphviz）。

**每一步的命令都会先写成脚本：** `<目录>/scripts/<编号>_<rule>[.<通配符>].sh`，然后再执行；日志写入 `<目录>/logs/`，软件版本写入 `<目录>/logs/versions/`。排查问题时，激活脚本头部注明的环境，直接 `bash` 该脚本即可单独重跑这一步。

## 快速开始

```bash
# 1. 配置
cp config/config.example.yaml config/config.yaml        # 英文注释版：config/config.example.en.yaml（参数相同）
#    必填：genome、rna.short / rna.long、homologs、qc.busco_lineages、qc.omark_db
#    功能注释：function.interproscan_dir、function.eggnog_data_dir
#    样本表：从 config/*.example.tsv 复制为 config/rna_short.tsv、rna_long.tsv、homologs.tsv 再填写

# 2. 外部代码
#   BRAKER4 是 git 子模块（external/BRAKER4，固定在 commit 6a5ee16），克隆本仓库时一并获取：
#   git clone --recursive <本流程仓库> annotation   （已克隆但缺少 BRAKER4 时：git submodule update --init）
#   （可选）git clone https://github.com/hyphaltip/EVidenceModeler_rust external/EVidenceModeler_rust

# 3. 运行
snakemake --configfile config/config.yaml --profile profiles/local -n                     # dry-run
snakemake --configfile config/config.yaml --profile profiles/local                        # 单机；按 workflow/envs 建 conda 环境（~/annotation_conda_envs）
snakemake --configfile config/config.yaml --profile profiles/slurm                        # SLURM（另设 braker4.executor: slurm）

# 只跑某一部分
snakemake --configfile config/config.yaml --profile profiles/local data_qc         # 00_Data 评估
snakemake --configfile config/config.yaml --profile profiles/local --until evm     # 到 EVM 为止
```

软件清单与版本见 `workflow/envs/README.md`，各环境的版本均已锁定。

### 在独立的分析目录中运行

建议把流程代码和每次分析分开：流程代码只放一份，每个物种 / 每次分析各用一个目录，放配置、样本表和结果。
用 `-d` 把 Snakemake 的工作目录设为分析目录，配置里的相对路径、`results/` 和 `.snakemake/` 都在分析目录内，不会写进代码目录：

```
my_species/
├── config.yaml          从 config/config.example.yaml 复制后修改
├── rna_short.tsv  rna_long.tsv  homologs.tsv
└── results/             运行后生成（00_Data ... 07_Report）
```

```bash
snakemake -s /path/to/annotation/workflow/Snakefile -d my_species --configfile my_species/config.yaml \
          --profile /path/to/annotation/profiles/local -n      # 去掉 -n 即正式运行
```

## 根据数据自动做出的选择

| 数据 | BRAKER4 | EviAnn 标签 | EVM 转录本证据 | GeMoMa RNA 证据 | 修正步骤的 RNA 支持 |
|---|---|---|---|---|---|
| 短读 + PacBio (+ONT) | dual | bam_mix | 长读 StringTie | 短读 BAM | 长读内含子链 |
| 短读 + ONT | ETP（不用 ONT） | bam_mix | 长读 StringTie | 短读 BAM | 长读内含子链 |
| 只有短读 | ETP | bam | 短读 StringTie | 短读 BAM | 短读：剪接 reads（B）+ StringTie 内含子链（A）+ 单外显子（C） |
| 只有 PacBio | IsoSeq | bam_isoseq | 长读 StringTie | 无（可选长读） | 长读内含子链 |
| 只有 ONT | EP | bam_isoseq（暂定） | 长读 StringTie | 无（可选长读） | 长读内含子链 |
| 没有 RNA（需 `rna.allow_none: true`） | EP | 不运行 | 无 | 无 | 无：修正只做 R5 证据检查；过滤以蛋白证据代替 RNA 证据；无 UTR / 可变剪接 / lncRNA |

这些选择的逻辑在 `workflow/lib/common.py`，测试在 `tests/test_common.py`。

## 关键规则

**EVM 权重调优（`05_EVM/tune`，配置 `evm.tune`）**
- 在一条序列（默认最长者）上，对各证据来源的候选权重做网格搜索，默认 16 组加配置中的基线。
- 用原始证据给每组结果打分：P 为多外显子模型的内含子链被长读 / 短读（与修正步骤同一标准）/ miniprot 任一完全支持的比例，R 为证据位点被覆盖的比例，F = 2PR/(P+R)。
- 默认用 F 最高的权重跑全基因组 EVM。

**修正（`05_EVM/refine`，`workflow/scripts/refine.py`）**
- 按 CDS 内含子链比较 EVM 和 EviAnn 模型。
- 长读支持以"全长 reads"计：Iso-Seq HQ 转录本按表头的 `full_length_coverage` 加权，原始 reads 按 1 计（`refine.lr_weight_by_coverage`）。
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
- CDS 短于 `min_cds_aa` 的模型单独剔除。
- TEsorter 结果只记录，不作为删除条件。
- 每个基因的判断写入 `filter_decisions.tsv`。

## 测试

```bash
python -m pytest tests                                                     # 单元测试
snakemake -s workflow/Snakefile --configfile tests/config.test.yaml -n     # 假数据 dry-run
```
