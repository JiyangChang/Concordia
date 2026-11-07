#!/usr/bin/env python3
"""Replace internal sequence IDs (column 1 of GFF3) by the original IDs from seqid_map.tsv."""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    m = {}
    with open(a.map) as fh:
        next(fh)
        for line in fh:
            new, orig = line.rstrip("\n").split("\t")
            m[new] = orig
    with open(a.gff) as fh, open(a.out, "w") as out:
        for line in fh:
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) > 1 and p[1] in m:
                    p[1] = m[p[1]]
                out.write(" ".join(p) + "\n")
                continue
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue
            cols = line.split("\t", 1)
            if cols[0] not in m:
                sys.exit(f"ERROR: sequence {cols[0]} not in seqid map")
            out.write(m[cols[0]] + "\t" + cols[1])


if __name__ == "__main__":
    main()
