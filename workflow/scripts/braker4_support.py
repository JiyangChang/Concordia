#!/usr/bin/env python3
"""Summarise BRAKER4 gene_support.tsv: for every evidence column, how many models are supported.

The file layout is not documented by BRAKER4, so every column after the first (the model ID) is
summarised: numeric columns -> % of models with value > 0 and mean; text columns -> % 'yes'-like.
An empty or missing file gives a one-line note instead.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from filter_final import _truthy, read_gene_support  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--support", required=True)
    ap.add_argument("--mode", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open(a.out, "w") as fh:
        fh.write("column\tmodels\tsupported\tsupported_%\tmean\n")
        if not os.path.exists(a.support) or os.path.getsize(a.support) == 0:
            fh.write(f"# BRAKER4 mode {a.mode}; gene_support.tsv missing or empty\n")
            return
        header, body = read_gene_support(a.support)
        fh.write(f"# BRAKER4 mode {a.mode}; {len(body)} models\n")
        for i, col in enumerate(header):
            if col in ("gene_id", "transcript_id", "chrom", "strand"):
                continue
            vals = [r[i] for r in body if i < len(r) and r[i] not in ("", "NA")]
            nums = []
            for v in vals:
                try:
                    nums.append(float(v.rstrip("%")))
                except ValueError:
                    pass
            sup = sum(_truthy(v) for v in vals)
            mean = f"{sum(nums) / len(nums):.3f}" if nums and len(nums) == len(vals) else "-"
            fh.write(f"{col}\t{len(vals)}\t{sup}\t{100 * sup / len(vals) if vals else 0:.1f}\t{mean}\n")


if __name__ == "__main__":
    main()
