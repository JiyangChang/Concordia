#!/usr/bin/env python3
"""Combine the per-chunk InterProScan results into one TSV and one GFF3.

  merge_interproscan.py --tsv chunk_001.tsv ... --gff chunk_001.gff3 ... --out-tsv all.tsv --out-gff all.gff3

TSV: concatenated (InterProScan TSV has no header). GFF3: one header (##gff-version, ##feature-ontology,
##interproscan-version from the first chunk that has them), then all ##sequence-region lines and features, then a
single ##FASTA section with the sequences of all chunks. Empty chunk files (no proteins) are allowed.
"""
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", nargs="*", default=[])
    ap.add_argument("--gff", nargs="*", default=[])
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-gff", required=True)
    a = ap.parse_args()

    n = 0
    with open(a.out_tsv, "w") as out:
        for p in a.tsv:
            for line in open(p):
                if line.strip():
                    out.write(line if line.endswith("\n") else line + "\n")
                    n += 1

    header, regions, features, fasta = [], [], [], []
    for p in a.gff:
        in_fasta = False
        for line in open(p):
            line = line.rstrip("\n")
            if not line:
                continue
            if in_fasta:
                fasta.append(line)
            elif line == "##FASTA":
                in_fasta = True
            elif line.startswith("##sequence-region"):
                regions.append(line)
            elif line.startswith("##"):
                key = line.split()[0]
                if key not in {h.split()[0] for h in header}:
                    header.append(line)
            else:
                features.append(line)
    if not any(h.startswith("##gff-version") for h in header):
        header.insert(0, "##gff-version 3")
    with open(a.out_gff, "w") as out:
        out.write("\n".join(header + regions + features) + "\n")
        if fasta:
            out.write("##FASTA\n" + "\n".join(fasta) + "\n")
    print(f"{len(a.tsv)} chunks: {n} TSV rows, {len(regions)} proteins, {len(features)} GFF3 features")


if __name__ == "__main__":
    main()
