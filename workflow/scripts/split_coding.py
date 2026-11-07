#!/usr/bin/env python3
"""Split an annotation into protein-coding genes and non-coding (lncRNA) genes."""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import parse_attrs  # noqa: E402
from gff_models import load_models, is_coding  # noqa: E402


def is_pseudo(g):
    """EviAnn marks processed pseudogenes with gene_biotype=processed_pseudogene;pseudo=true."""
    at = parse_attrs(g["row"][8])
    return at.get("pseudo", "").lower() == "true" or "pseudogene" in at.get("gene_biotype", "").lower()


def write_gene(fh, g):
    fh.write("\t".join(g["row"]) + "\n")
    for t in g["tx"].values():
        fh.write("\t".join(t["row"]) + "\n")
        for f in t["feats"]:
            fh.write("\t".join(f) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--coding", required=True)
    ap.add_argument("--noncoding", required=True)
    ap.add_argument("--pseudogenes", default="", help="pseudogenes are written here (never to --coding)")
    a = ap.parse_args()
    genes, _ = load_models(a.gff)
    nc = c = ps = 0
    fp = open(a.pseudogenes, "w") if a.pseudogenes else None
    if fp:
        fp.write("##gff-version 3\n")
    with open(a.coding, "w") as fc, open(a.noncoding, "w") as fn:
        fc.write("##gff-version 3\n")
        fn.write("##gff-version 3\n")
        for g in genes.values():
            if is_pseudo(g):
                ps += 1
                if fp:
                    write_gene(fp, g)
            elif is_coding(g):
                write_gene(fc, g)
                c += 1
            else:
                write_gene(fn, g)
                nc += 1
    if fp:
        fp.close()
    print(f"coding genes: {c}\tnon-coding genes: {nc}\tpseudogenes: {ps}")


if __name__ == "__main__":
    main()
