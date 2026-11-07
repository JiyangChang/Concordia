#!/usr/bin/env python3
"""Write BRAKER4 samples.csv and config.ini for one genome."""
import argparse
import configparser
import csv
import os

COLUMNS = ["sample_name", "genome", "genome_masked", "protein_fasta", "bam_files", "fastq_r1", "fastq_r2",
           "sra_ids", "varus_genus", "varus_species", "isoseq_bam", "isoseq_fastq", "busco_lineage",
           "reference_gtf"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--masked", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--bam", default="")
    ap.add_argument("--isoseq-bam", default="")
    ap.add_argument("--busco-lineage", required=True)
    ap.add_argument("--min-contig", default="10000")
    ap.add_argument("--skip-busco", default="1")
    ap.add_argument("--run-omark", default="0")
    ap.add_argument("--use-minisplice", default="0")
    ap.add_argument("--threads", default="8")
    ap.add_argument("--mem-mb", default="64000")
    ap.add_argument("--runtime", default="4320")
    ap.add_argument("--template", default="", help="BRAKER4 repo config.ini used as the base (keeps its [paths]/[containers])")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-ini", required=True)
    a = ap.parse_args()

    row = {c: "" for c in COLUMNS}
    row.update(sample_name=a.sample, genome=a.genome, genome_masked=a.masked, protein_fasta=a.proteins,
               bam_files=a.bam, isoseq_bam=a.isoseq_bam, busco_lineage=a.busco_lineage)
    with open(a.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerow(row)
    ini = configparser.ConfigParser()
    ini.optionxform = str  # keep key case
    if a.template and os.path.exists(a.template):
        ini.read(a.template)
    for sec in ("PARAMS", "SLURM_ARGS"):
        if not ini.has_section(sec):
            ini.add_section(sec)
    ini["PARAMS"].update(min_contig=a.min_contig, skip_busco=a.skip_busco, run_omark=a.run_omark,
                         use_minisplice=a.use_minisplice)
    ini["SLURM_ARGS"].update(cpus_per_task=a.threads, mem_of_node=a.mem_mb, max_runtime=a.runtime)
    # Point BRAKER4 at our samples.csv if its [paths] section has a samples entry
    if ini.has_section("paths"):
        for k in list(ini["paths"].keys()):
            if "sample" in k.lower():
                ini["paths"][k] = os.path.abspath(a.out_csv)
    with open(a.out_ini, "w") as fh:
        ini.write(fh)

if __name__ == "__main__":
    main()
