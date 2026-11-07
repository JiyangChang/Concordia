#!/usr/bin/env python3
"""Keep one transcript per gene: the one with the longest CDS (ties: longest exon sum, then the first).

Lines are copied unchanged except that the gene row is shrunk to the extent of the kept transcript.
Features of dropped transcripts are removed; a feature with several parents keeps only the kept one.
Comment and directive lines are kept.

  select_longest.py --gff all_isoforms.gff3 --out longest.gff3 [--report report.tsv]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import parse_attrs  # noqa: E402
from gff_models import TX_TYPES, is_gene_type  # noqa: E402


def set_attr_parent(field, parents):
    parts = []
    for p in field.rstrip(";").split(";"):
        if p.startswith("Parent="):
            p = "Parent=" + ",".join(parents)
        parts.append(p)
    return ";".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report")
    a = ap.parse_args()

    with open(a.gff) as fh:
        lines = fh.read().splitlines()
    rows = [l.split("\t") if l and not l.startswith("#") else None for l in lines]

    tx_gene, cds_len, exon_len, tx_span, order = {}, {}, {}, {}, {}
    for c in rows:
        if c is None or len(c) < 9:
            continue
        at = parse_attrs(c[8])
        if c[2] in TX_TYPES:
            tid = at["ID"]
            tx_gene[tid] = at.get("Parent", "").split(",")[0]
            tx_span[tid] = (int(c[3]), int(c[4]))
            order[tid] = len(order)
    for c in rows:
        if c is None or len(c) < 9 or c[2] not in ("CDS", "exon"):
            continue
        n = int(c[4]) - int(c[3]) + 1
        for p in parse_attrs(c[8]).get("Parent", "").split(","):
            d = cds_len if c[2] == "CDS" else exon_len
            d[p] = d.get(p, 0) + n

    best = {}
    for tid, gid in tx_gene.items():
        key = (cds_len.get(tid, 0), exon_len.get(tid, 0), -order[tid])
        if gid not in best or key > best[gid][0]:
            best[gid] = (key, tid)
    keep = {tid for _, tid in best.values()}
    n_tx_genes = {}
    for tid, gid in tx_gene.items():
        n_tx_genes[gid] = n_tx_genes.get(gid, 0) + 1

    out = []
    for line, c in zip(lines, rows):
        if c is None or len(c) < 9:
            out.append(line)
            continue
        at = parse_attrs(c[8])
        if is_gene_type(c[2]) and at.get("ID") in best:
            s, e = tx_span[best[at["ID"]][1]]
            c = c[:3] + [str(s), str(e)] + c[5:]
            out.append("\t".join(c))
        elif c[2] in TX_TYPES:
            if at["ID"] in keep:
                out.append(line)
        elif "Parent" in at:
            parents = at["Parent"].split(",")
            if not any(p in tx_gene for p in parents):     # child of a gene etc.: keep
                out.append(line)
                continue
            kept = [p for p in parents if p in keep]
            if kept:
                out.append(line if kept == parents else "\t".join(c[:8] + [set_attr_parent(c[8], kept)]))
        else:
            out.append(line)
    with open(a.out, "w") as fh:
        fh.write("\n".join(out) + "\n")

    multi = sum(1 for n in n_tx_genes.values() if n > 1)
    rep = (f"genes\t{len(best)}\ntranscripts_in_input\t{len(tx_gene)}\ntranscripts_kept\t{len(keep)}\n"
           f"genes_with_several_isoforms\t{multi}\n")
    if a.report:
        with open(a.report, "w") as fh:
            fh.write("item\tvalue\n" + rep)
    sys.stdout.write(rep)


if __name__ == "__main__":
    main()
