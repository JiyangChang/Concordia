#!/usr/bin/env python3
"""Cut EVM inputs down to one sequence for weight tuning (standard library only).

--region auto  -> the longest sequence; otherwise an ORIGINAL sequence name (e.g. Chr1), mapped to
the internal ID through seqid_map.tsv. Writes <outdir>/genome.fa, every given GFF3 restricted to
that sequence (same file names), and region.tsv (internal id, original name, length).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_fasta, write_fasta  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--seqid-map", required=True)
    ap.add_argument("--region", default="auto")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("gffs", nargs="+")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    to_new, to_orig = {}, {}
    with open(a.seqid_map) as fh:
        next(fh)
        for line in fh:
            new, orig = line.rstrip("\n").split("\t")
            to_new[orig], to_orig[new] = new, orig
    seqs = {h.split()[0]: s for h, s in read_fasta(a.genome)}
    if a.region == "auto":
        region = max(seqs, key=lambda k: len(seqs[k]))
    else:
        region = to_new.get(a.region, a.region if a.region in seqs else None)
        if region is None:
            sys.exit(f"ERROR: tuning region {a.region} not found in the genome")
    with open(os.path.join(a.outdir, "genome.fa"), "w") as fh:
        write_fasta(fh, region, seqs[region])
    for g in a.gffs:
        n = 0
        with open(g) as src, open(os.path.join(a.outdir, os.path.basename(g)), "w") as dst:
            for line in src:
                if line.startswith("#") or line.split("\t", 1)[0] == region:
                    dst.write(line)
                    n += not line.startswith("#")
        print(f"{os.path.basename(g)}: {n} features on {region}")
    with open(os.path.join(a.outdir, "region.tsv"), "w") as fh:
        fh.write("internal_id\toriginal_id\tlength\n")
        fh.write(f"{region}\t{to_orig.get(region, region)}\t{len(seqs[region])}\n")


if __name__ == "__main__":
    main()
