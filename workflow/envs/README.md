# 软件环境

版本用 conda（libmamba 求解器，conda-forge + bioconda）求解并锁定直接依赖的确切版本；Python 等基础依赖交给求解器决定。
2026-11-03 首次锁定；2026-11-07 把 evm、assembly_qc、readqc、report 并入 core 后重新求解。

默认建立 7 个环境；merqury、evm-rust 只在启用对应功能时建立。

| 环境 | 用途 | 主要版本 |
|---|---|---|
| core | 输入整理、比对与组装、GeMoMa、miniprot、EVM、端粒 / 污染筛查、长读质控、报告 | HISAT2 2.2.3, minimap2 2.31, samtools 1.24, StringTie 3.0.3, GeMoMa 1.9, miniprot 0.18, gffread 0.12.9, gffcompare 0.12.10, fastp 1.4.0, chopper 0.15.0, seqkit 2.14.0, EVidenceModeler 2.1.0, tidk 0.3.1, NCBI FCS-GX 0.5.5, NanoPlot 1.46.2, matplotlib 3.10.9 |
| repeat | 重复序列、LAI、TEsorter | EDTA 2.3.0, TEsorter 1.5.1, pandas 2.3.3 |
| eviann | RNA-based 注释 | EviAnn 2.0.7 |
| agat | GFF 处理与统计 | AGAT 1.7.0, gffread 0.12.9 |
| qc | 蛋白 / 基因组完整度 | BUSCO 6.1.0, MultiQC 1.35 |
| omark | 蛋白组一致性 | OMArk 0.5.0, OMAmer 2.1.2 |
| function | 功能注释 | eggNOG-mapper 2.1.12, openjdk 11（InterProScan 用官方安装包） |
| merqury（可选） | k-mer 评估 | Merqury 1.4.1, meryl 1.4.2 |
| evm-rust（可选） | Rust 版 EVM（源码编译） | rust 1.99.1 |

为什么不再合并（用 conda 求解器和实际安装验证过）：
- **eviann**：EviAnn 的包自带 samtools、stringtie、gffread、gffcompare、blastp、trmap 等可执行文件和与 EVM 同名的 `bin/PerlLib`，与 core 装在一起会互相覆盖（实测 core 中的 samtools 变成 1.15.1、stringtie 变成 3.0.0）。求解器检查不出这类文件冲突。
- **repeat**：EDTA 内部的 TIR-Learner 3.0.7 不兼容 pandas 3，必须用 pandas 2.3.3；环境本身也最大。
- **agat**：AGAT 1.7.0 依赖旧的 R 4.1，与 core 的库冲突，无法求解。
- **qc**：BUSCO 6.1 依赖多（AUGUSTUS、MetaEuk 等），与 NanoPlot、OMArk 冲突。
- **omark / function**：两者可以合并求解，但 eggNOG-mapper 的 diamond 会从 2.2.8 降到 2.1.10、pandas 也随之改变，收益（1–2 GB）不值得重新核对结果。

BRAKER4 没有单独的环境：嵌套运行 BRAKER4 时直接使用启动流程的 Snakemake（`annotation` 环境：Snakemake 8.18.2 + snakemake-executor-plugin-slurm 2.6.0），BRAKER4 的工具都在它自带的容器中。
