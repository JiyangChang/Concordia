# 更新记录

[English](CHANGELOG.md) | 中文（中英文两版内容相同，修改时请同步）

## v0.1.0（2026-11-20）

首个发布版本。

### 流程

- **输入：** 基因组组装；RNA-seq 可以是短读和 / 或长读（PacBio Iso-Seq / HiFi、ONT cDNA / dRNA），也可以完全没有
  （`rna.allow_none: true`）；同源物种的基因组和注释；OrthoDB 蛋白。
- **00_Data：** 输入检查（序列名、同源注释只保留蛋白编码基因）、组装质量（连续性、BUSCO 基因组模式、LAI、端粒；可选
  Merqury 和 FCS-GX）、reads 质控。
- **01_Repeat：** EDTA 构建 TE 库并注释，RepeatMasker 软屏蔽，重复序列统计。
- **02_RNA：** HISAT2 / minimap2 比对、StringTie 组装、EviAnn 注释（每个样本一行输入，可用 `tissue` / `sample_group`
  列合并重复样本），输出编码基因和 lncRNA。
- **03_Ab_initio：** BRAKER4（git 子模块，嵌套 Snakemake，Apptainer 容器），模式按数据自动选择：dual / ETP / IsoSeq / EP。
- **04_Homology：** 每个同源物种单独运行 GeMoMa（再用 GAF 合并），以及 miniprot 蛋白比对。
- **05_EVM：** EVidenceModeler 整合（可做权重调优，默认关闭）；与 EviAnn 比较修正（R1–R6 规则：UTR、可变剪接、补入 /
  替换 / 拆分 / 融合基因，RNA 支持按长读或短读判断）；TE 与长度过滤（有 RNA 证据的短基因保留）；按染色体命名；导出主结果
  （每个基因一个转录本）和全部异构体。
- **06_Function：** InterProScan（分块运行，合并为一个 TSV / GFF3）和 eggNOG-mapper；GO、KEGG、COG 统计。
- **每个注释来源的质控：** BUSCO、OMArk、AGAT；可选与参考注释比较。
- **07_Report：** 中英文两份独立的 HTML 报告，每一步都有表和图；各注释来源对比（BUSCO 和 OMArk 图）、软件版本，以及每个
  作业的运行时间和资源用量，并给出 `resources` 的调整建议。
- 每一步都把命令写成脚本（`<步骤>/scripts/`），并记录日志、软件版本和 benchmark。
- 下载目录中缺少的 BUSCO 数据集自动下载；按自身中间文件断点续跑的工具（EDTA、EviAnn、EVidenceModeler）只在命令和输入
  都没变时续跑。
- 可在单机或 SLURM 上运行（profile 在 `profiles/`，完整示例见 `docs/example_slurm.zh.md`）；启动时检查 SLURM 的 `cores`
  设置和 `default-resources` 的值。

### 运行要求

- Snakemake 8.18.2 与 snakemake-executor-plugin-slurm 2.6.0（BRAKER4 测试过的版本），`setuptools<81`。
- BRAKER4 子模块固定在 commit 6a5ee16；Apptainer（或 Singularity 3.x）。
- conda 环境按 `workflow/envs/*.yaml` 建立（默认 7 个）。
- 运行节点的 CPU 需支持 AVX2（bioconda 的 StringTie 构建需要）。
- 数据库：见 `scripts/download_databases.sh` 和 `docs/deployment_guide.zh.md`。

### 已测试

- 辣木（*Moringa oleifera*）13 号染色体（13.2 Mb），工作站上从头运行，分别使用四种 RNA 数据：短读 + 长读、仅短读、
  仅长读、无 RNA-seq。
- 辣木全基因组（236 Mb），SLURM 集群：20,394 个蛋白编码基因、772 个 lncRNA，BUSCO 98.3%（eudicotyledons_odb12），
  92.5% 的基因有功能注释。

### 已知限制

- EDTA / EviAnn / EVidenceModeler 的断点续跑检查在只改线程数时也会判定为"已变化"，工具会从头开始。
- `gemoma.ere_long_reads`（用长读作为 GeMoMa 的 RNA 证据）为实验性选项，默认关闭。
- 只有 ONT 数据时以 `bam_isoseq` 交给 EviAnn（暂定方案）；ONT reads 不进入 BRAKER4。
- EviAnn 只输出至少在 2 个 RNA 样本中检测到的 lncRNA。
- "运行时间与资源"一节不统计 BRAKER4 容器内的进程；请以 SLURM 的 `sacct` 为准。
