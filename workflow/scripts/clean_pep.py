#!/usr/bin/env python3
"""Clean a protein FASTA for BUSCO / OMArk / InterProScan / eggNOG-mapper: strip terminal stop
symbols ('*' or '.'), replace internal stops with 'X' (InterProScan rejects them) and drop empty
records. The number of proteins with internal stops is reported on stderr."""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_fasta, write_fasta  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("inp")
ap.add_argument("out")
a = ap.parse_args()
n_internal = 0
with open(a.out, "w") as fh:
    for h, s in read_fasta(a.inp):
        s = s.rstrip("*.")
        if "*" in s or "." in s:
            n_internal += 1
            s = s.replace("*", "X").replace(".", "X")
        if s:
            write_fasta(fh, h.split()[0], s)
if n_internal:
    print(f"clean_pep: {n_internal} protein(s) with internal stop codons; internal stops written as X",
          file=sys.stderr)
