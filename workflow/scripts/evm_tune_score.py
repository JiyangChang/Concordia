#!/usr/bin/env python3
"""Score EVM runs of a weight grid on the tuning region and pick the best weights.

For each combination (EVM GFF3 of the region):
  precision P: fraction of multi-exon transcripts whose CDS intron chain is confirmed by ANY
    available evidence -- long reads (same rule as refine: read introns inside the CDS span == CDS
    introns; weighted by full-length reads), or, only without long reads, the short-read support of
    refine.py (1 intron: junction reads, absolute and relative to competing introns; >= 2 introns:
    that plus an identical StringTie intron chain), or an identical miniprot alignment chain.
    Requiring any (not each) evidence type avoids penalising genes that only have RNA support
    (e.g. lineage-specific genes without homologs).
    Per-evidence rates (lr_chain, sr_chain, prot_chain) are reported for information.
  recall R: mean of
    lr_loci    fraction of spliced long reads overlapping an EVM CDS (either strand)
    prot_loci  fraction of miniprot alignments with Identity >= --min-identity overlapping an EVM
               CDS on the same strand
  F = 2PR / (P + R).
The best combination has the highest F (ties: higher P, then lower combination number). Its
weights file is copied to --best. Note: miniprot and transcripts are also EVM inputs, so the score
slightly favours those evidence types.
"""
import argparse
import os
import shutil
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import parse_attrs, read_gff  # noqa: E402
from refine import Index, LongReadSupport, ShortReadSupport, genes_from  # noqa: E402


def miniprot(path, region, min_identity):
    """(set of (seqid, strand, CDS intron chain), [(seqid, strand, s, e)] of good alignments)."""
    ident, segs = {}, defaultdict(list)
    for c in read_gff(path):
        if c[0] != region:
            continue
        at = parse_attrs(c[8])
        if c[2] == "mRNA":
            try:
                ident[at["ID"]] = float(at.get("Identity", 0))
            except ValueError:
                ident[at["ID"]] = 0.0
        elif c[2] == "CDS":
            segs[(at.get("Parent", ""), c[0], c[6])].append((int(c[3]), int(c[4])))
    chains, good = set(), []
    for (pid, seqid, strand), sg in segs.items():
        sg = sorted(sg)
        chain = tuple((x[1] + 1, y[0] - 1) for x, y in zip(sg, sg[1:]))
        if chain:
            chains.add((seqid, strand, chain))
        if ident.get(pid, 0.0) >= min_identity:
            good.append((seqid, strand, sg[0][0], sg[-1][1]))
    return chains, good


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combo", action="append", required=True, help="ID=evm.gff3")
    ap.add_argument("--weights-dir", required=True, help="contains <ID>/weights.txt")
    ap.add_argument("--region-file", required=True)
    ap.add_argument("--miniprot", required=True)
    ap.add_argument("--lr-chains", default="")
    ap.add_argument("--sr-junctions", default="")
    ap.add_argument("--min-lr-support", type=int, default=2)
    ap.add_argument("--min-sr-junction", type=int, default=3)
    ap.add_argument("--sr-min-junction-ratio", type=float, default=0.1)
    ap.add_argument("--sr-stringtie-min-cov", type=float, default=3.0)
    ap.add_argument("--stringtie", default="", help="StringTie GTF (short-read criterion A)")
    ap.add_argument("--min-identity", type=float, default=0.8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--best", required=True)
    a = ap.parse_args()

    with open(a.region_file) as fh:
        next(fh)
        region, orig, length = fh.readline().rstrip("\n").split("\t")
    lr = LongReadSupport(a.lr_chains, 0.5) if a.lr_chains and os.path.exists(a.lr_chains) else None
    sr = (ShortReadSupport(a.sr_junctions, a.stringtie, a.min_sr_junction, a.sr_min_junction_ratio,
                           a.sr_stringtie_min_cov, 0.5)
          if lr is None and a.sr_junctions and os.path.exists(a.sr_junctions) else None)
    prot_chains, prot_good = miniprot(a.miniprot, region, a.min_identity)
    spliced_reads = [r for r in (lr.reads if lr else []) if r[0] == region and r[3]]

    rows = []
    for item in a.combo:
        cid, path = item.split("=", 1)
        genes = genes_from(path, "EVM") if os.path.exists(path) else []
        txs = [t for g in genes for t in g.txs]
        multi = [t for t in txs if len(t.cds) > 1]
        ok = {}
        if lr is not None:
            ok["lr_chain"] = [lr.support(t.seqid, t.cds) >= a.min_lr_support for t in multi]
        if sr is not None:
            ok["sr_chain"] = [sr.support(t.seqid, t.cds, t.strand) > 0 for t in multi]
        ok["prot_chain"] = [(t.seqid, t.strand, t.cds_introns) in prot_chains for t in multi]
        m = {k: (sum(v) / len(v) if v else 0.0) for k, v in ok.items()}
        confirmed = [any(v[i] for v in ok.values()) for i in range(len(multi))]
        cds_any = Index([(t.seqid, s, e) for t in txs for s, e in t.cds])
        cds_str = Index([((t.seqid, t.strand), s, e) for t in txs for s, e in t.cds])
        r = {}
        if spliced_reads:
            r["lr_loci"] = sum(cds_any.covered(sq, s, e) > 0 for sq, s, e, _ in spliced_reads) / len(spliced_reads)
        if prot_good:
            r["prot_loci"] = sum(cds_str.covered((sq, st), s, e) > 0 for sq, st, s, e in prot_good) / len(prot_good)
        P = sum(confirmed) / len(confirmed) if confirmed else 0.0
        R = sum(r.values()) / len(r) if r else 0.0
        F = 2 * P * R / (P + R) if P + R else 0.0
        w = open(os.path.join(a.weights_dir, cid, "weights.txt")).read().strip().splitlines()
        wstr = ", ".join(f"{x.split(chr(9))[1]}={x.split(chr(9))[2]}" for x in w)
        rows.append((cid, wstr, len(genes), len(multi), m, r, P, R, F))

    rows.sort(key=lambda x: (-x[8], -x[6], x[0]))
    mk = sorted({k for x in rows for k in x[4]})
    rk = sorted({k for x in rows for k in x[5]})
    with open(a.out, "w") as fh:
        fh.write(f"# region {orig} ({region}, {length} bp); P = multi-exon chains confirmed by any of "
                 f"{', '.join(mk)}; R = mean({', '.join(rk)}), "
                 f"F = 2PR/(P+R); c00 = configured evm.weights\n")
        fh.write("\t".join(["rank", "combo", "weights", "genes", "multi_exon_tx"] + mk + rk + ["P", "R", "F"]) + "\n")
        for i, (cid, wstr, ng, nm, m, r, P, R, F) in enumerate(rows, 1):
            fh.write("\t".join(map(str, [i, cid, wstr, ng, nm] + [f"{m.get(k, 0):.4f}" for k in mk] +
                                   [f"{r.get(k, 0):.4f}" for k in rk] + [f"{P:.4f}", f"{R:.4f}", f"{F:.4f}"])) + "\n")
    shutil.copy(os.path.join(a.weights_dir, rows[0][0], "weights.txt"), a.best)
    print(f"best: {rows[0][0]} ({rows[0][1]}) F={rows[0][8]:.4f}")


if __name__ == "__main__":
    main()
