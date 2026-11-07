#!/usr/bin/env python3
"""Collect BUSCO / OMArk / AGAT results of every annotation source into one table."""
import argparse
import glob
import os
import re

BUSCO_RE = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)")


def busco(path):
    for f in glob.glob(os.path.join(path, "short_summary*.txt")):
        m = BUSCO_RE.search(open(f).read())
        if m:
            return dict(zip(["C", "S", "D", "F", "M", "n"], m.groups()))
    return {}


def omark(path):
    out = {}
    for f in glob.glob(os.path.join(path, "*.sum")):
        for line in open(f):
            if line.startswith("S:"):
                out["omark_completeness"] = line.strip()
            elif line.startswith("A:"):
                out["omark_consistency"] = line.strip()
    return out


def agat(path):
    out = {}
    f = os.path.join(path, "stats.txt")
    if os.path.exists(f):
        for line in open(f):
            m = re.match(r"\s*Number of gene\s+(\d+)", line)
            if m and "genes" not in out:
                out["genes"] = m.group(1)
            m = re.match(r"\s*Number of mrna\s+(\d+)", line)
            if m and "mrnas" not in out:
                out["mrnas"] = m.group(1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--lineages", nargs="+", required=True)
    ap.add_argument("entries", nargs="+", help="name=qc_dir")
    a = ap.parse_args()
    cols = ["source", "lineage", "C", "S", "D", "F", "M", "n", "omark_completeness", "omark_consistency", "genes", "mrnas"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in a.entries:
            name, d = e.split("=", 1)
            base = {"source": name, **omark(os.path.join(d, "omark")), **agat(os.path.join(d, "agat"))}
            for lin in a.lineages:
                row = {**base, "lineage": lin, **busco(os.path.join(d, "busco", lin))}
                fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


if __name__ == "__main__":
    main()
