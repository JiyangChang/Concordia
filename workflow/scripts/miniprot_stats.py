#!/usr/bin/env python3
"""Per-species miniprot alignment statistics.

--species name=proteins.faa,miniprot.gff3 (repeatable)
Columns: species, input proteins, proteins with >= 1 alignment (%), alignments (mRNA records),
mean identity, alignments with identity >= 0.9, complete alignments (start M and stop codon found:
miniprot marks frameshifts/stops in the attributes, so 'complete' = no StopCodon/Frameshift issue).
"""
import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_gff, parse_attrs  # noqa: E402


def count_fasta(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", action="append", required=True, help="name=proteins.faa,miniprot.gff3")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cols = ["species", "input_proteins", "aligned_proteins", "aligned_%", "alignments",
            "mean_identity", "identity>=0.9", "with_frameshift_or_stop"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for item in a.species:
            name, files = item.split("=", 1)
            faa, gff = files.split(",", 1)
            n_in = count_fasta(faa)
            targets, idents, n_aln, n_bad = set(), [], 0, 0
            for c in read_gff(gff):
                if c[2] != "mRNA":
                    continue
                n_aln += 1
                at = parse_attrs(c[8])
                targets.add(at.get("Target", "").split(" ")[0])
                try:
                    idents.append(float(at.get("Identity", "nan")))
                except ValueError:
                    pass
                if int(at.get("Frameshift", 0) or 0) > 0 or int(at.get("StopCodon", 0) or 0) > 0:
                    n_bad += 1
            idents = [x for x in idents if x == x]
            row = [name, n_in, len(targets), f"{100 * len(targets) / n_in:.1f}" if n_in else "-", n_aln,
                   f"{statistics.mean(idents):.3f}" if idents else "-",
                   sum(x >= 0.9 for x in idents), n_bad]
            fh.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
