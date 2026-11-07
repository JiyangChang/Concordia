# Eukaryotic genome structural and functional annotation workflow (Snakemake)

English | [中文](README.md)

Given a target genome, RNA-seq (short reads and/or PacBio / ONT long reads) and the genomes and annotations of
several homologous species, the workflow runs the following steps and ends with an illustrated statistics report:
- assembly quality assessment
- repeat annotation and soft-masking
- RNA alignment
- three gene prediction routes (EviAnn / BRAKER4 / GeMoMa)
- EVM integration
- refinement against EviAnn
- filtering
- functional annotation

- Deployment guide: `docs/deployment_guide.en.md` (Chinese version `docs/deployment_guide.md`)
- Example for a SLURM cluster: `docs/example_slurm.en.md` (Chinese version `docs/example_slurm.md`)
- The Chinese and English documents (README, deployment guide, SLURM example, configuration template, profiles) have the same
  content; when editing, update both versions.

## Overview

| Directory | Content | Main tools |
|---|---|---|
| `00_Data` | Input collection (optionally a BRAKER4 protein subset by OrthoDB level, `orthodb_subset`; unchanged inputs are symbolic links, not copies: genome, single homolog protein / CDS files, the protein database when only OrthoDB is used); sequences renamed when a name exceeds 13 characters; homolog annotations reduced to protein-coding genes; assembly QC (contiguity, BUSCO, telomeres, LAI (after EDTA), optional Merqury / FCS-GX); raw RNA data QC | seqkit, BUSCO, tidk, fastp, NanoPlot |
| `01_Repeat` | TE library and soft-masking; TE composition, Kimura landscape, LTR insertion times | EDTA, RepeatMasker |
| `02_RNA` | Short- and long-read alignment, merged BAMs, StringTie assembly, long-read intron chains; EviAnn annotation (RNA data only; SwissProt for CDS decisions, homolog proteins if not configured) | HISAT2, minimap2, samtools, StringTie2, EviAnn |
| `03_Ab_initio` | BRAKER4 (nested Snakemake; uses PacBio long reads only) | BRAKER4 |
| `04_Homology` | Homology projection and protein alignment: GeMoMa and miniprot per homologous species (each with BUSCO / OMArk), then merged (GeMoMa GAF; miniprot concatenated) for EVM | GeMoMa, miniprot |
| `05_EVM` | EVM integration → refinement (R1–R6) → TE / length filter → chromosome-based names (e.g. `Mol01G000010`) → original sequence names restored → export; BUSCO / OMArk / AGAT of every source | EVidenceModeler 2.1 (Rust version optional), AGAT, TEsorter |
| `06_Function` | Domains, GO, KEGG, COG | InterProScan, eggNOG-mapper |
| `07_Report` | `annotation_report.html`: tables and figures for every step | Python |

**Workflow diagrams:**
- `docs/workflow_main.png` / `.svg`: main workflow (with RNA-seq) grouped by the 00–07 directories, without QC and statistics steps; an arrow means "reads the output of", and every direct dependency is drawn.
- `docs/workflow_main_norna.png` / `.svg`: main workflow in no-RNA mode (`rna.allow_none: true`).
- `docs/workflow_rulegraph.png` / `.svg`: all steps (103 with the example configuration; QC and statistics steps in light grey).
- Regenerate after changing rules: `bash scripts/make_workflow_diagrams.sh` (needs snakemake and graphviz; uses the example configuration `tests/config.test.yaml`, no databases needed).

**The command of every step is first written to a script:** `<directory>/scripts/<number>_<rule>[.<wildcards>].sh`, which
is then executed; logs go to `<directory>/logs/`, software versions to `<directory>/logs/versions/`. To debug a step,
activate the environment named in the script header and run the script with `bash`.

## Quick start

```bash
# 1. Configuration
cp config/config.example.en.yaml config/config.yaml     # Chinese comments: config/config.example.yaml (same parameters)
#    required: genome, rna.short / rna.long, homologs, qc.busco_lineages, qc.omark_db
#    functional annotation: function.interproscan_dir, function.eggnog_data_dir
#    sample tables: copy config/*.example.tsv to config/rna_short.tsv, rna_long.tsv, homologs.tsv and fill them in

# 2. External code
#   BRAKER4 is a git submodule (external/BRAKER4, pinned to commit 6a5ee16), fetched together with this repository:
#   git clone --recursive <this workflow's repository> annotation   (already cloned without BRAKER4: git submodule update --init)
#   (optional) git clone https://github.com/hyphaltip/EVidenceModeler_rust external/EVidenceModeler_rust

# 3. Run
snakemake --configfile config/config.yaml --profile profiles/local -n                     # dry-run
snakemake --configfile config/config.yaml --profile profiles/local                        # single machine; conda environments from workflow/envs (~/annotation_conda_envs)
snakemake --configfile config/config.yaml --profile profiles/slurm                        # SLURM (also set braker4.executor: slurm)

# Run only a part
snakemake --configfile config/config.yaml --profile profiles/local data_qc         # 00_Data assessment
snakemake --configfile config/config.yaml --profile profiles/local --until evm     # up to EVM
```

Software list and versions: `workflow/envs/README.md`; the versions of all environments are pinned.

### Running in a separate analysis directory

Keep the workflow code and each analysis apart: one copy of the code, and one directory per species / analysis with
its configuration, sample tables and results. With `-d`, Snakemake's working directory is the analysis directory:
relative paths in the configuration, `results/` and `.snakemake/` all stay there and never touch the code directory:

```
my_species/
├── config.yaml          copied from config/config.example.en.yaml and edited
├── rna_short.tsv  rna_long.tsv  homologs.tsv
└── results/             created by the run (00_Data ... 07_Report)
```

```bash
snakemake -s /path/to/annotation/workflow/Snakefile -d my_species --configfile my_species/config.yaml \
          --profile /path/to/annotation/profiles/local -n      # drop -n for the real run
```

## Choices made from the data

| Data | BRAKER4 | EviAnn tag | EVM transcript evidence | GeMoMa RNA evidence | RNA support in refinement |
|---|---|---|---|---|---|
| short + PacBio (+ONT) | dual | bam_mix | long-read StringTie | short-read BAM | long-read intron chains |
| short + ONT | ETP (ONT not used) | bam_mix | long-read StringTie | short-read BAM | long-read intron chains |
| short reads only | ETP | bam | short-read StringTie | short-read BAM | short reads: junction reads (B) + StringTie intron chain (A) + single exon (C) |
| PacBio only | IsoSeq | bam_isoseq | long-read StringTie | none (long reads optional) | long-read intron chains |
| ONT only | EP | bam_isoseq (provisional) | long-read StringTie | none (long reads optional) | long-read intron chains |
| no RNA (needs `rna.allow_none: true`) | EP | not run | none | none | none: refinement only checks R5 evidence; the filter uses protein evidence instead of RNA; no UTRs / alternative isoforms / lncRNAs |

The logic is in `workflow/lib/common.py`, tested in `tests/test_common.py`.

## Key rules

**EVM weight tuning (`05_EVM/tune`, configuration `evm.tune`)**
- Grid search over candidate weights of the evidence sources on one sequence (the longest by default): 16
  combinations by default plus the configured baseline.
- Each result is scored with the raw evidence: P = fraction of multi-exon models whose intron chain is fully
  supported by long reads / short reads (same criteria as refinement) / miniprot (any), R = fraction of evidence
  loci covered, F = 2PR/(P+R).
- By default the weights with the highest F are used for the genome-wide EVM run.

**Refinement (`05_EVM/refine`, `workflow/scripts/refine.py`)**
- EVM and EviAnn models are compared by their CDS intron chains.
- Long-read support counts "full-length reads": Iso-Seq HQ transcripts are weighted by `full_length_coverage` from
  the header, raw reads count 1 (`refine.lr_weight_by_coverage`).
- Rules: R1 imports UTRs and isoforms; R2 adds missing genes; R3 replaces a conflicting model with the one supported
  by long reads; R4 handles splits and fusions; R5 checks the evidence of EVM-only genes; R6 lncRNAs are output
  separately.
- The decision for every locus is written to `refine_decisions.tsv`. Loci that need manual review are listed
  separately in the report.

**Naming (`05_EVM/final`, `workflow/scripts/rename_genes.py`, configuration `naming`)**
- Coding genes `<prefix><chromosome>G<number>`, e.g. `Mol01G000010`; lncRNAs use `L`, e.g. `Mol01L000010`.
- Numbered by position within each chromosome, step 10.
- The chromosome number is taken from the original sequence name (Chr1 / chr01 / Chromosome_3 / 1); unmatched
  scaffolds share `00`; codes can be set by hand, e.g. `ChrM: Mt`.
- Transcripts are `gene.1`, `.2` (`.1` = longest CDS).
- The main result `<prefix>.gff3` (and cds/pep/cdna) keeps one transcript per gene (longest CDS; on ties the longest
  transcript); all isoforms are kept in `<prefix>.all_isoforms.*`. Statistics and assessments of the final result
  are based on the main result.
- 5′/3′ UTR features are added; the old -> new ID table is `id_map.tsv`.

**Filtering (`05_EVM/final`)**
- A gene is removed only if its CDS overlaps TEs by more than `max_te_overlap` (default 0.5) **and** it has no RNA
  evidence.
- RNA evidence, any of: StringTie, from EviAnn, BRAKER4 gene_support.
- Models with a CDS shorter than `min_cds_aa` are removed separately.
- TEsorter hits are only recorded, never a removal criterion.
- The decision for every gene is written to `filter_decisions.tsv`.

## Tests

```bash
python -m pytest tests                                                     # unit tests
snakemake -s workflow/Snakefile --configfile tests/config.test.yaml -n     # dry-run on mock data
```
