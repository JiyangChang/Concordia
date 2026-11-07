#!/usr/bin/env python3
"""Split a FASTA into exactly N chunk files (round-robin); chunks may be empty."""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_fasta, write_fasta  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("inp")
ap.add_argument("outs", nargs="+")
a = ap.parse_args()
fhs = [open(p, "w") for p in a.outs]
for i, (h, s) in enumerate(read_fasta(a.inp)):
    write_fasta(fhs[i % len(fhs)], h, s)
for fh in fhs:
    fh.close()
