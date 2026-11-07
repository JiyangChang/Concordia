#!/usr/bin/env python3
"""Agreement of annotation sets with a reference annotation, at the CDS level.

The reference is restricted to sequences present in the (internal) genome; its sequence IDs are
translated through seqid_map.tsv (original -> internal). CDS-level comparison is used because
many reference annotations (e.g. AUGUSTUS-based) carry no UTRs.

Per set (name=gff3):
  gene_recall      fraction of reference genes whose CDS overlaps a predicted gene's CDS (same strand)
  gene_precision   fraction of predicted genes whose CDS overlaps a reference gene's CDS
  intron_recall / intron_precision   CDS introns (seqid, strand, start, end) shared
  tx_recall        fraction of reference transcripts whose complete CDS structure is predicted exactly
  tx_precision     fraction of predicted transcripts with an exact reference CDS structure
A reference that is itself a prediction measures agreement, not accuracy.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff_models import load_models  # noqa: E402
from refine import Index  # noqa: E402


def load(path, rename=None, keep=None):
    """[(gene_id, seqid, strand, [cds tuples per transcript])] for coding genes."""
    out = []
    models, _ = load_models(path)
    for gid, g in models.items():
        seqid = g["row"][0]
        if rename is not None:
            seqid = rename.get(seqid)
        if seqid is None or (keep is not None and seqid not in keep):
            continue
        txs = []
        for t in g["tx"].values():
            cds = tuple(sorted((int(f[3]), int(f[4])) for f in t["feats"] if f[2] == "CDS"))
            if cds:
                txs.append(cds)
        if txs:
            out.append((gid, seqid, g["row"][6], txs))
    return out


def summarise(ref, pred):
    def index(genes):
        return Index([((sq, st), s, e) for _, sq, st, txs in genes for cds in txs for s, e in cds])

    def introns(genes):
        return {(sq, st, a[1] + 1, b[0] - 1) for _, sq, st, txs in genes for cds in txs for a, b in zip(cds, cds[1:])}

    def chains(genes):
        return {(sq, st, cds) for _, sq, st, txs in genes for cds in txs}

    def hits(genes, ix):
        return sum(any(ix.covered((sq, st), s, e) > 0 for cds in txs for s, e in cds) for _, sq, st, txs in genes)

    ri, pi = index(ref), index(pred)
    rint, pint = introns(ref), introns(pred)
    rch, pch = chains(ref), chains(pred)
    frac = lambda a, b: f"{100 * a / b:.1f}" if b else "-"
    return {
        "ref_genes": len(ref), "pred_genes": len(pred),
        "gene_recall_%": frac(hits(ref, pi), len(ref)),
        "gene_precision_%": frac(hits(pred, ri), len(pred)),
        "intron_recall_%": frac(len(rint & pint), len(rint)),
        "intron_precision_%": frac(len(rint & pint), len(pint)),
        "tx_recall_%": frac(len(rch & pch), len(rch)),
        "tx_precision_%": frac(len(rch & pch), len(pch)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--seqid-map", required=True)
    ap.add_argument("--set", action="append", required=True, help="name=gff3 (internal sequence IDs)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    orig2new = {}
    with open(a.seqid_map) as fh:
        next(fh)
        for line in fh:
            new, orig = line.rstrip("\n").split("\t")
            orig2new[orig] = new
    ref = load(a.reference, rename=orig2new)
    cols = ["set", "ref_genes", "pred_genes", "gene_recall_%", "gene_precision_%", "intron_recall_%",
            "intron_precision_%", "tx_recall_%", "tx_precision_%"]
    with open(a.out, "w") as fh:
        fh.write(f"# reference: {os.path.basename(a.reference)}; {len(ref)} coding genes on the analysed sequences; "
                 "CDS-level comparison (agreement, not accuracy, if the reference is itself a prediction)\n")
        fh.write("\t".join(cols) + "\n")
        for item in a.set:
            name, path = item.split("=", 1)
            pred = load(path) if os.path.exists(path) else []
            row = {"set": name, **summarise(ref, pred)}
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    print(open(a.out).read())


if __name__ == "__main__":
    main()
