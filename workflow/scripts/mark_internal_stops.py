#!/usr/bin/env python3
"""Mark transcripts whose CDS translation contains internal stop codons (they are kept).

Adds internal_stop=yes to the transcript rows (and to their gene) of --gff, rewritten in place,
using the protein FASTA translated from it (gffread -y: stops are '.' or '*').
--report writes the marked transcripts.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_fasta, parse_attrs  # noqa: E402
from gff_models import TX_TYPES, is_gene_type  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--pep", required=True)
    ap.add_argument("--report")
    a = ap.parse_args()
    bad = {h.split()[0] for h, s in read_fasta(a.pep) if "*" in s.rstrip("*.") or "." in s.rstrip("*.")}
    lines = open(a.gff).read().splitlines()
    genes = set()
    for l in lines:
        c = l.split("\t")
        if len(c) >= 9 and c[2] in TX_TYPES:
            at = parse_attrs(c[8])
            if at.get("ID") in bad:
                genes.add(at.get("Parent", "").split(",")[0])
    out = []
    for l in lines:
        c = l.split("\t")
        if len(c) >= 9 and not l.startswith("#"):
            at = parse_attrs(c[8])
            if ((c[2] in TX_TYPES and at.get("ID") in bad) or (is_gene_type(c[2]) and at.get("ID") in genes)) \
                    and "internal_stop=" not in c[8]:
                c[8] = c[8].rstrip(";") + ";internal_stop=yes"
                l = "\t".join(c)
        out.append(l)
    tmp = a.gff + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\n".join(out) + "\n")
    os.replace(tmp, a.gff)
    if a.report:
        with open(a.report, "w") as fh:
            fh.write("transcript\tnote\n" + "".join(f"{t}\tinternal stop codon(s); kept, marked internal_stop=yes\n"
                                                     for t in sorted(bad)))
    print(f"transcripts with internal stop codons: {len(bad)} (genes: {len(genes)})")


if __name__ == "__main__":
    main()
