#!/usr/bin/env python3
"""Summarise RNA-seq alignments (per sample and per merged BAM) for 02_RNA/qc and the report.

Per sample (from samtools stats SN lines, spliced count and, for short reads, the HISAT2 summary):
  reads, mapped, mapping_%, unique_%, multi_% (HISAT2 only), spliced_% of mapped reads,
  error_rate (mismatches / mapped bases, needs NM tags), mean read length.
Per merged BAM (from samtools flagstat): total, mapped, mapping_%, and which steps use it.
"""
import argparse
import os
import re


def sn(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                out[p[0].rstrip(":")] = p[1]
    return out


def hisat2(path):
    """overall %, unique % and multi % (concordant for PE, overall for SE)."""
    if not os.path.exists(path):
        return {}
    t = open(path).read()
    d = {}
    m = re.search(r"([\d.]+)% overall alignment rate", t)
    if m:
        d["overall"] = m.group(1)
    m = re.search(r"\(([\d.]+)%\) aligned concordantly exactly 1 time", t) or re.search(r"\(([\d.]+)%\) aligned exactly 1 time", t)
    if m:
        d["unique"] = m.group(1)
    m = re.search(r"\(([\d.]+)%\) aligned concordantly >1 times", t) or re.search(r"\(([\d.]+)%\) aligned >1 times", t)
    if m:
        d["multi"] = m.group(1)
    return d


def flagstat(path):
    d = {}
    if os.path.exists(path):
        for line in open(path):
            m = re.match(r"(\d+) \+ \d+ in total", line)
            if m:
                d["total"] = int(m.group(1))
            m = re.match(r"(\d+) \+ \d+ primary mapped \(([\d.]+)%", line)
            if m:
                d["primary_mapped"], d["primary_pct"] = int(m.group(1)), m.group(2)
            m = re.match(r"(\d+) \+ \d+ primary$", line.strip()) or re.match(r"(\d+) \+ \d+ primary\s*$", line)
            if m:
                d["primary"] = int(m.group(1))
            m = re.match(r"(\d+) \+ \d+ mapped \(([\d.]+)%", line)
            if m:
                d["mapped"], d["mapped_pct"] = int(m.group(1)), m.group(2)
    return d


def pct(a, b):
    return f"{100 * a / b:.2f}" if b else "0"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--samples", nargs="*", default=[], help="kind:sample:platform")
    ap.add_argument("--merged", nargs="*", default=[], help="group=use1;use2")
    ap.add_argument("--out-samples", required=True)
    ap.add_argument("--out-merged", required=True)
    a = ap.parse_args()

    cols = ["sample", "type", "platform", "reads", "mapped", "mapping_%", "unique_%", "multi_%",
            "spliced_%", "error_rate", "mean_length"]
    with open(a.out_samples, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in a.samples:
            kind, name, platform = e.split(":")
            s = sn(os.path.join(a.qc_dir, kind, f"{name}.samtools_stats.txt"))
            total = int(s.get("raw total sequences", 0))
            mapped = int(s.get("reads mapped", 0))
            sp_file = os.path.join(a.qc_dir, kind, f"{name}.spliced.txt")
            spliced = int(open(sp_file).read().strip() or 0) if os.path.exists(sp_file) else 0
            h = hisat2(os.path.join(a.qc_dir, f"{name}.hisat2.txt")) if kind == "short" else {}
            row = [name, kind, platform, total, mapped, pct(mapped, total), h.get("unique", "-"), h.get("multi", "-"),
                   pct(spliced, mapped), s.get("error rate", "-"), s.get("average length", "-")]
            fh.write("\t".join(map(str, row)) + "\n")

    with open(a.out_merged, "w") as fh:
        fh.write("bam\ttotal_records\tprimary_mapped\tprimary_mapped_%\tused_by\n")
        for e in a.merged:
            g, _, uses = e.partition("=")
            f = flagstat(os.path.join(a.qc_dir, f"{g}.flagstat.txt"))
            pm = f.get("primary_mapped", f.get("mapped", 0))
            pp = f.get("primary_pct", f.get("mapped_pct", "0"))
            fh.write(f"{g}.bam\t{f.get('total', 0)}\t{pm}\t{pp}\t{uses.replace(';', '; ') or '-'}\n")


if __name__ == "__main__":
    main()
