#!/usr/bin/env python3
"""Gene-model statistics for one or more GFF3 sets (standard library only).

Usage: annotation_stats.py --set coding=a.gff3 --set lncRNA=b.gff3 --out stats.tsv --out-dist dist.tsv

stats.tsv  metric x set table:
  genes, transcripts, transcripts_per_gene, multi_isoform_genes(_%), exons_per_transcript,
  single_exon_transcripts(_%), gene/transcript/CDS length (mean, median), total exon bp,
  coding transcripts, transcripts with 5' UTR / 3' UTR / both (_%).
dist.tsv   long table (set, kind, value) for plots: isoforms_per_gene, exons_per_transcript,
  transcript_length, cds_aa.
UTRs are taken from explicit five_prime_UTR / three_prime_UTR features when present, otherwise
derived from exons extending beyond the CDS (strand-aware).
"""
import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff_models import load_models  # noqa: E402

UTR5 = {"five_prime_UTR", "5UTR", "five_prime_utr"}
UTR3 = {"three_prime_UTR", "3UTR", "three_prime_utr"}


def tx_stats(t, strand):
    feats = t["feats"]
    exons = sorted((int(f[3]), int(f[4])) for f in feats if f[2] == "exon")
    cds = sorted((int(f[3]), int(f[4])) for f in feats if f[2] == "CDS")
    if not exons:
        exons = cds
    tx_len = sum(e - s + 1 for s, e in exons)
    cds_bp = sum(e - s + 1 for s, e in cds)
    types = {f[2] for f in feats}
    has5 = bool(types & UTR5)
    has3 = bool(types & UTR3)
    if cds and not (types & (UTR5 | UTR3)):
        lo, hi = cds[0][0], cds[-1][1]
        left = any(s < lo for s, _ in exons)
        right = any(e > hi for _, e in exons)
        has5, has3 = (left, right) if strand != "-" else (right, left)
    return {"exons": len(exons), "len": tx_len, "cds_aa": cds_bp // 3, "coding": bool(cds),
            "utr5": has5, "utr3": has3}


def mean(x):
    return f"{statistics.mean(x):.1f}" if x else "0"


def median(x):
    return f"{statistics.median(x):.0f}" if x else "0"


def pct(a, b):
    return f"{100 * a / b:.1f}" if b else "0"


def summarise(path):
    genes, _ = load_models(path)
    iso, gene_len, txs = [], [], []
    for g in genes.values():
        if not g["tx"]:
            continue
        iso.append(len(g["tx"]))
        gene_len.append(int(g["row"][4]) - int(g["row"][3]) + 1)
        txs += [tx_stats(t, g["row"][6]) for t in g["tx"].values()]
    n_g, n_t = len(iso), len(txs)
    coding = [t for t in txs if t["coding"]]
    single = sum(t["exons"] == 1 for t in txs)
    multi = sum(i > 1 for i in iso)
    u5 = sum(t["utr5"] for t in coding)
    u3 = sum(t["utr3"] for t in coding)
    both = sum(t["utr5"] and t["utr3"] for t in coding)
    m = [
        ("genes", n_g), ("transcripts", n_t),
        ("transcripts_per_gene", f"{n_t / n_g:.2f}" if n_g else "0"),
        ("multi_isoform_genes", multi), ("multi_isoform_genes_%", pct(multi, n_g)),
        ("exons_per_transcript_mean", mean([t["exons"] for t in txs])),
        ("single_exon_transcripts", single), ("single_exon_transcripts_%", pct(single, n_t)),
        ("gene_length_mean", mean(gene_len)), ("gene_length_median", median(gene_len)),
        ("transcript_length_mean", mean([t["len"] for t in txs])),
        ("transcript_length_median", median([t["len"] for t in txs])),
        ("coding_transcripts", len(coding)),
        ("cds_length_aa_mean", mean([t["cds_aa"] for t in coding])),
        ("cds_length_aa_median", median([t["cds_aa"] for t in coding])),
        ("coding_tx_with_5UTR_%", pct(u5, len(coding))), ("coding_tx_with_3UTR_%", pct(u3, len(coding))),
        ("coding_tx_with_both_UTR_%", pct(both, len(coding))),
    ]
    dist = ([("isoforms_per_gene", i) for i in iso] + [("exons_per_transcript", t["exons"]) for t in txs]
            + [("transcript_length", t["len"]) for t in txs] + [("cds_aa", t["cds_aa"]) for t in coding])
    return m, dist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", action="append", required=True, help="name=path.gff3")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-dist", required=True)
    a = ap.parse_args()
    names, results = [], []
    with open(a.out_dist, "w") as fd:
        fd.write("set\tkind\tvalue\n")
        for item in a.set:
            name, path = item.split("=", 1)
            m, dist = summarise(path)
            names.append(name)
            results.append(dict(m))
            for kind, v in dist:
                fd.write(f"{name}\t{kind}\t{v}\n")
    metrics = list(results[0].keys())  # same order for every set
    with open(a.out, "w") as fh:
        fh.write("metric\t" + "\t".join(names) + "\n")
        for k in metrics:
            fh.write(k + "\t" + "\t".join(str(r.get(k, "-")) for r in results) + "\n")


if __name__ == "__main__":
    main()
