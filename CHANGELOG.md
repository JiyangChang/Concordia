# Changelog

English | [中文](CHANGELOG.zh.md) (both versions have the same content - keep them in sync when editing)

## v0.1.0 (2026-11-20)

First release.

### Workflow

- **Inputs:** a genome assembly; RNA-seq as short reads and/or long reads (PacBio Iso-Seq / HiFi, ONT cDNA / dRNA) or
  none at all (`rna.allow_none: true`); genomes and annotations of homologous species; OrthoDB proteins.
- **00_Data:** input checks (sequence names, homolog annotations cleaned to protein-coding genes), assembly quality
  (contiguity, BUSCO genome mode, LAI, telomeres; optional Merqury and FCS-GX), read QC.
- **01_Repeat:** EDTA TE library and annotation, RepeatMasker soft-masking, repeat statistics.
- **02_RNA:** HISAT2 / minimap2 alignment, StringTie assemblies, EviAnn annotation with one input line per sample
  (optional `tissue` / `sample_group` column to combine replicates), coding genes and lncRNAs.
- **03_Ab_initio:** BRAKER4 (git submodule, nested Snakemake, Apptainer containers), mode chosen from the data:
  dual / ETP / IsoSeq / EP.
- **04_Homology:** GeMoMa per homologous species (combined with GAF) and miniprot protein alignments.
- **05_EVM:** EVidenceModeler integration (weight tuning available, off by default), refinement against EviAnn
  (rules R1-R6: UTRs, isoforms, added / replaced / split / fused genes, RNA support from long or short reads), TE and
  length filter (short genes with RNA evidence kept), chromosome-based gene names, export of the main result (one
  transcript per gene) and of all isoforms.
- **06_Function:** InterProScan (chunked, merged into one TSV / GFF3) and eggNOG-mapper; GO, KEGG, COG statistics.
- **QC of every annotation source:** BUSCO, OMArk, AGAT; optional comparison with a reference annotation.
- **07_Report:** self-contained HTML report in English and Chinese with tables and figures for every step, a
  comparison of all annotation sources (BUSCO and OMArk figures), software versions, and the runtime and resources of
  every job with suggestions for the `resources` settings.
- Every step writes its command to a script (`<step>/scripts/`), logs, software versions and benchmark files.
- BUSCO datasets missing from the download directory are downloaded automatically; tools that resume from their own
  intermediate files (EDTA, EviAnn, EVidenceModeler) only resume when command and inputs are unchanged.
- Runs locally or on SLURM (profiles in `profiles/`, worked example in `docs/example_slurm.md`); start-up checks for
  the SLURM `cores` setting and the `default-resources` values.

### Requirements

- Snakemake 8.18.2 with snakemake-executor-plugin-slurm 2.6.0 (the versions BRAKER4 is tested with), `setuptools<81`.
- BRAKER4 submodule at commit 6a5ee16; Apptainer (or Singularity 3.x).
- conda environments are built from `workflow/envs/*.yaml` (7 by default).
- Nodes whose CPU supports AVX2 (required by the bioconda StringTie build).
- Databases: see `scripts/download_databases.sh` and `docs/deployment_guide.md`.

### Tested

- *Moringa oleifera* chromosome 13 (13.2 Mb) on a workstation, from scratch, with each kind of RNA data: short and
  long reads, short reads only, long reads only, and no RNA-seq.
- *Moringa oleifera* whole genome (236 Mb) on a SLURM cluster: 20,394 protein-coding genes, 772 lncRNAs, BUSCO
  98.3% (eudicotyledons_odb12), 92.5% of genes with functional annotation.

### Known limitations

- The resume check of EDTA / EviAnn / EVidenceModeler also changes when only the number of threads changes, so the
  tool then starts again from the beginning.
- `gemoma.ere_long_reads` (long reads as GeMoMa RNA evidence) is experimental and off by default.
- ONT-only data are given to EviAnn as `bam_isoseq` (provisional); ONT reads are not passed to BRAKER4.
- EviAnn reports lncRNAs only when they are found in at least 2 RNA samples.
- The runtime and resources section does not count the processes inside the BRAKER4 containers; use SLURM's `sacct`.
