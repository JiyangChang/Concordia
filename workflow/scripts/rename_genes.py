#!/usr/bin/env python3
"""Chromosome-based gene names for the final annotation (standard library only).

Coding genes : <prefix><chrom><G><number>   e.g. Mol01G000010
lncRNA genes : <prefix><chrom><L><number>   e.g. Mol01L000010
  chrom  : chromosome number taken from the ORIGINAL sequence name with --chrom-regex (Chr1,
           chr01, Chromosome_3, 1 -> 01), zero-padded to max(2, digits of the largest number);
           --code SEQID=CODE overrides it (e.g. ChrM=Mt); sequences that match neither are
           unplaced and share --unplaced-code (default 00), numbered in genome order.
  number : per chromosome code, genes ordered by position, numbered step, 2*step, ...
           (default 10, 20, ...) with --digits digits, leaving room for later insertions.
--scheme sequential gives <prefix><G><number> numbered over the whole genome instead.
Transcripts are <gene>.1, <gene>.2, ... (.1 = longest CDS; for lncRNA the longest transcript).
Writes gene / mRNA (lnc_RNA) / exon / CDS / five_prime_UTR / three_prime_UTR with recomputed CDS
phases, keeps the refine_rule / eviann_support / refine_support gene attributes, records the old
ID as source_id=, and writes an old -> new ID table.
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import parse_attrs  # noqa: E402
from gff_models import load_models  # noqa: E402

KEEP_ATTRS = ("refine_rule", "eviann_support", "refine_support")


def seq_order(map_path):
    """internal seqid -> (index in genome order, original seqid)."""
    order = {}
    with open(map_path) as fh:
        next(fh)
        for i, line in enumerate(fh):
            new, orig = line.rstrip("\n").split("\t")
            order[new] = (i, orig)
    return order


def chrom_codes(order, regex, overrides, unplaced):
    rx = re.compile(regex, re.I)
    num = {}
    for new, (_, orig) in order.items():
        if orig in overrides:
            continue
        m = rx.match(orig)
        if m:
            num[new] = int(m.group(1))
    width = max(2, len(str(max(num.values())))) if num else 2
    seen = {}
    for new, n in num.items():
        if n in seen:
            sys.exit(f"ERROR: sequences {order[seen[n]][1]} and {order[new][1]} both map to chromosome {n}; "
                     "use naming.chrom_codes to disambiguate")
        seen[n] = new
    codes = {}
    for new, (_, orig) in order.items():
        codes[new] = overrides.get(orig) or (f"{num[new]:0{width}d}" if new in num else unplaced)
    return codes


def structure(t):
    feats = t["feats"]
    cds = sorted((int(f[3]), int(f[4])) for f in feats if f[2] == "CDS")
    exons = sorted((int(f[3]), int(f[4])) for f in feats if f[2] == "exon") or list(cds)
    return exons, cds


def write_models(fh, genes, new_ids, kind, step_rows, order):
    """kind = 'mRNA' (coding) or 'lnc_RNA'. Genes are written in genome order (sequence, position)."""
    pos = lambda g: (order[genes[g]["row"][0]][0], int(genes[g]["row"][3]), int(genes[g]["row"][4]))
    for gid in sorted(genes, key=pos):
        g = genes[gid]
        new = new_ids[gid]
        row = g["row"]
        seqid, strand, src = row[0], row[6], row[1]
        at = parse_attrs(row[8])
        txs = []
        for tid, t in g["tx"].items():
            exons, cds = structure(t)
            if exons:
                txs.append((tid, exons, cds))
        if kind == "mRNA":
            txs.sort(key=lambda x: (-sum(e - s + 1 for s, e in x[2]), x[1][0][0]))
        else:
            txs.sort(key=lambda x: (-sum(e - s + 1 for s, e in x[1]), x[1][0][0]))
        gs = min(x[1][0][0] for x in txs)
        ge = max(x[1][-1][1] for x in txs)
        extra = "".join(f";{k}={at[k]}" for k in KEEP_ATTRS if k in at)
        biotype = "protein_coding" if kind == "mRNA" else "lncRNA"
        fh.write(f"{seqid}\t{src}\tgene\t{gs}\t{ge}\t.\t{strand}\t.\t"
                 f"ID={new};Name={new};biotype={biotype};source_id={gid}{extra}\n")
        step_rows.append(("gene", gid, new, seqid, gs, ge))
        for k, (tid, exons, cds) in enumerate(txs, 1):
            tnew = f"{new}.{k}"
            step_rows.append((kind, tid, tnew, seqid, exons[0][0], exons[-1][1]))
            fh.write(f"{seqid}\t{src}\t{kind}\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t"
                     f"ID={tnew};Parent={new};Name={tnew};source_id={tid}\n")
            tx_order = exons if strand != "-" else list(reversed(exons))
            for i, (a, b) in enumerate(tx_order, 1):
                fh.write(f"{seqid}\t{src}\texon\t{a}\t{b}\t.\t{strand}\t.\tID={tnew}.exon{i};Parent={tnew}\n")
            if not cds:
                continue
            lo, hi = cds[0][0], cds[-1][1]
            left = [(a, min(b, lo - 1)) for a, b in exons if a < lo]
            right = [(max(a, hi + 1), b) for a, b in exons if b > hi]
            utr5, utr3 = (left, right) if strand != "-" else (right, left)
            for typ, segs, tagname in (("five_prime_UTR", utr5, "utr5p"), ("three_prime_UTR", utr3, "utr3p")):
                for i, (a, b) in enumerate(segs, 1):
                    fh.write(f"{seqid}\t{src}\t{typ}\t{a}\t{b}\t.\t{strand}\t.\tID={tnew}.{tagname}{i};Parent={tnew}\n")
            done, phase = 0, {}
            for a, b in (cds if strand != "-" else list(reversed(cds))):
                phase[(a, b)] = (3 - done % 3) % 3
                done += b - a + 1
            for a, b in cds:
                fh.write(f"{seqid}\t{src}\tCDS\t{a}\t{b}\t.\t{strand}\t{phase[(a, b)]}\tID={tnew}.cds;Parent={tnew}\n")


def assign(genes, order, codes, prefix, tag, scheme, step, digits):
    def key(gid):
        r = genes[gid]["row"]
        return (order.get(r[0], (10 ** 9, r[0]))[0], int(r[3]), int(r[4]), gid)

    missing = sorted({genes[g]["row"][0] for g in genes if genes[g]["row"][0] not in order})
    if missing:
        sys.exit(f"ERROR: sequences not in seqid map: {missing[:5]}")
    new_ids, counters = {}, {}
    for gid in sorted(genes, key=key):
        code = "" if scheme == "sequential" else codes[genes[gid]["row"][0]]
        counters[code] = counters.get(code, 0) + step
        if counters[code] >= 10 ** digits:
            sys.exit(f"ERROR: more genes on {code or 'genome'} than --digits {digits} allows with step {step}")
        new_ids[gid] = f"{prefix}{code}{tag}{counters[code]:0{digits}d}"
    return new_ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coding", required=True)
    ap.add_argument("--lncrna", default="")
    ap.add_argument("--seqid-map", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--scheme", choices=["chromosome", "sequential"], default="chromosome")
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--digits", type=int, default=6)
    ap.add_argument("--chrom-regex", default=r"^(?:chr(?:omosome)?)?[_-]?0*(\d+)$")
    ap.add_argument("--unplaced-code", default="00")
    ap.add_argument("--code", action="append", default=[], help="SEQID=CODE (original sequence name)")
    ap.add_argument("--coding-tag", default="G")
    ap.add_argument("--lncrna-tag", default="L")
    ap.add_argument("--out-coding", required=True)
    ap.add_argument("--out-lncrna", default="")
    ap.add_argument("--out-map", required=True)
    a = ap.parse_args()

    order = seq_order(a.seqid_map)
    overrides = dict(x.split("=", 1) for x in a.code)
    codes = chrom_codes(order, a.chrom_regex, overrides, a.unplaced_code)
    rows = []

    coding, _ = load_models(a.coding)
    coding = {g: v for g, v in coding.items() if v["tx"]}
    ids = assign(coding, order, codes, a.prefix, a.coding_tag, a.scheme, a.step, a.digits)
    with open(a.out_coding, "w") as fh:
        fh.write("##gff-version 3\n")
        write_models(fh, coding, ids, "mRNA", rows, order)

    n_lnc = 0
    if a.lncrna and a.out_lncrna:
        lnc, _ = load_models(a.lncrna) if os.path.getsize(a.lncrna) else ({}, [])
        lnc = {g: v for g, v in lnc.items() if v["tx"]}
        lids = assign(lnc, order, codes, a.prefix, a.lncrna_tag, a.scheme, a.step, a.digits)
        with open(a.out_lncrna, "w") as fh:
            fh.write("##gff-version 3\n")
            write_models(fh, lnc, lids, "lnc_RNA", rows, order)
        n_lnc = len(lnc)

    with open(a.out_map, "w") as fh:
        fh.write("feature\told_id\tnew_id\tseqid\tstart\tend\n")
        for kind, old, new, seqid, s, e in rows:
            fh.write(f"{kind}\t{old}\t{new}\t{order[seqid][1]}\t{s}\t{e}\n")
    used = sorted({codes[v["row"][0]] for v in coding.values()}) if a.scheme == "chromosome" else []
    print(f"coding genes: {len(coding)}; lncRNA genes: {n_lnc}; chromosome codes used: {', '.join(used) or '-'}")


if __name__ == "__main__":
    main()
