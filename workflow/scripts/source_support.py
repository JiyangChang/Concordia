#!/usr/bin/env python3
"""How do the transcripts of one annotation (EVM / final) relate to the input sources?

For every coding transcript of --query and every --source name=gff:
  exact   : the source has a transcript with the identical CDS chain (seqid, strand, all CDS
            segment coordinates)
  overlap : at least one CDS base overlaps a CDS of the source on the same strand
Outputs
  --out-support : source, transcripts, exact, exact_%, overlap, overlap_%  (+ a 'none' row:
                  transcripts overlapping no source at all)
  --out-combos  : combination of sources with an exact match -> number of transcripts
                  (the data behind an UpSet-style bar chart)
"""
import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff_models import load_models  # noqa: E402
from filter_final import IntervalIndex  # noqa: E402


def chains(path):
    """List of (key, cds segments) for coding transcripts; key = (seqid, strand, chain)."""
    genes, _ = load_models(path)
    out = []
    for g in genes.values():
        seqid, strand = g["row"][0], g["row"][6]
        for t in g["tx"].values():
            cds = tuple(sorted((int(f[3]), int(f[4])) for f in t["feats"] if f[2] == "CDS"))
            if cds:
                out.append(((seqid, strand, cds), cds))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True)
    ap.add_argument("--source", action="append", required=True, help="name=path.gff3")
    ap.add_argument("--out-support", required=True)
    ap.add_argument("--out-combos", required=True)
    a = ap.parse_args()

    query = chains(a.query)
    names, exact_sets, indexes = [], [], []
    for item in a.source:
        name, path = item.split("=", 1)
        src = chains(path) if os.path.exists(path) else []
        names.append(name)
        exact_sets.append({k for k, _ in src})
        indexes.append(IntervalIndex([((k[0], k[1]), s, e) for k, cds in src for s, e in cds]))

    n = len(query)
    exact = Counter()
    overlap = Counter()
    combos = Counter()
    none = 0
    for key, cds in query:
        hit_exact = [nm for nm, es in zip(names, exact_sets) if key in es]
        hit_ov = [nm for nm, ix in zip(names, indexes) if any(ix.overlaps((key[0], key[1]), s, e) for s, e in cds)]
        exact.update(hit_exact)
        overlap.update(hit_ov)
        combos["+".join(hit_exact) if hit_exact else "no exact match"] += 1
        none += not hit_ov

    pct = lambda x: f"{100 * x / n:.1f}" if n else "0"
    with open(a.out_support, "w") as fh:
        fh.write("source\ttranscripts\texact\texact_%\toverlap\toverlap_%\n")
        for nm in names:
            fh.write(f"{nm}\t{n}\t{exact[nm]}\t{pct(exact[nm])}\t{overlap[nm]}\t{pct(overlap[nm])}\n")
        fh.write(f"none (no source overlaps)\t{n}\t-\t-\t{none}\t{pct(none)}\n")
    with open(a.out_combos, "w") as fh:
        fh.write("combination\ttranscripts\tpercent\n")
        for combo, c in combos.most_common():
            fh.write(f"{combo}\t{c}\t{pct(c)}\n")


if __name__ == "__main__":
    main()
