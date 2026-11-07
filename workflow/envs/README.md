# 软件环境

版本已于 2026-11-03 用 mamba（conda-forge + bioconda）求解并锁定直接依赖的确切版本；Python 等基础依赖交给求解器决定。

| 环境 | 用途 | 主要版本 |
|---|---|---|
| core | 比对、组装、GeMoMa、miniprot、gff 工具 | HISAT2 2.2.3, minimap2 2.31, samtools 1.24, StringTie 3.0.3, GeMoMa 1.9, miniprot 0.18, gffread 0.12.9, fastp 1.4.0 |
| repeat | 重复序列 | EDTA 2.3.0, TEsorter 1.5.1, pandas 2.3.3（TIR-Learner 不兼容 pandas 3） |
| eviann | RNA-based 注释 | EviAnn 2.0.7 |
| braker4 | 嵌套运行 BRAKER4 | Snakemake 8.18.2 + snakemake-executor-plugin-slurm 2.6.0（工具在 BRAKER4 自带容器中） |
| evm | 整合 | EVidenceModeler 2.1.0 |
| evm-rust | 可选的 Rust 版 EVM（源码编译） | rust 1.99.1 |
| agat | GFF 处理与统计 | AGAT 1.7.0 |
| qc | 蛋白 / 基因组完整度 | BUSCO 6.1.0, MultiQC 1.35 |
| omark | 蛋白组一致性 | OMArk 0.5.0, OMAmer 2.1.2 |
| readqc | 长读质控 | NanoPlot 1.46.2 |
| assembly_qc | 端粒、污染 | tidk 0.3.1, NCBI FCS-GX 0.5.5 |
| merqury | k-mer 评估 | Merqury 1.4.1, meryl 1.4.2 |
| function | 功能注释 | eggNOG-mapper 2.1.12（InterProScan 用官方安装包或容器） |
| report | 统计报告 | matplotlib 3.11.2, pandas 3.1.0 |

拆分说明：BUSCO 6.1 与 NanoPlot、OMArk 0.5 的依赖冲突（OMArk 需要 Python 3.11，BUSCO 6.1 需要 3.13），因此分成 qc / readqc / omark 三个环境。
