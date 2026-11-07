#!/usr/bin/env python3
"""Convert evidence into EVidenceModeler input formats.

Subcommands
  genepred   gene models -> gene/mRNA/exon/CDS with a fixed source column and source-prefixed IDs
  miniprot   miniprot --gff output -> PROTEIN alignments (nucleotide_to_protein_match)
  transcripts  StringTie GTF -> TRANSCRIPT alignments (cDNA_match with Target)
  repeats    RepeatMasker GFF -> simple repeat features with unique IDs
  weights    write weights.txt from "CLASS<TAB>SOURCE<TAB>WEIGHT" lines
"""
import argparse
import os
import sys
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_gff, parse_attrs  # noqa: E402
from gff_models import load_models  # noqa: E402


def cmd_genepred(a):
    genes, _ = load_models(a.gff)
    src, pre = a.source, a.source + "."
    out = open(a.out, "w")
    n = 0
    for gid, g in genes.items():
        coding_tx = [(tid, t) for tid, t in g["tx"].items() if any(f[2] == "CDS" for f in t["feats"])]
        if not coding_tx:
            continue
        r = g["row"][:]
        r[1], r[2], r[8] = src, "gene", f"ID={pre}{gid}"
        out.write("\t".join(r) + "\n")
        for tid, t in coding_tx:
            tr = t["row"][:]
            tr[1], tr[2], tr[8] = src, "mRNA", f"ID={pre}{tid};Parent={pre}{gid}"
            out.write("\t".join(tr) + "\n")
            exons = sorted((f for f in t["feats"] if f[2] == "exon"), key=lambda f: int(f[3]))
            cds = sorted((f for f in t["feats"] if f[2] == "CDS"), key=lambda f: int(f[3]))
            if not exons:  # e.g. GeMoMa: derive exons from CDS
                exons = [f[:] for f in cds]
            for i, e in enumerate(exons, 1):
                e = e[:]
                e[1], e[2], e[7], e[8] = src, "exon", ".", f"ID={pre}{tid}.exon{i};Parent={pre}{tid}"
                out.write("\t".join(e) + "\n")
            for i, cf in enumerate(cds, 1):
                cf = cf[:]
                cf[1], cf[2], cf[8] = src, "CDS", f"ID={pre}{tid}.cds;Parent={pre}{tid}"
                out.write("\t".join(cf) + "\n")
            n += 1
    out.close()
    print(f"{src}: {n} coding transcripts written")


def cmd_miniprot(a):
    mrna_target, n = {}, 0
    rows = list(read_gff(a.gff))
    for c in rows:
        if c[2] == "mRNA":
            at = parse_attrs(c[8])
            mrna_target[at["ID"]] = at.get("Target", "").split()[0] if at.get("Target") else at["ID"]
    with open(a.out, "w") as out:
        for c in rows:
            if c[2] != "CDS":
                continue
            at = parse_attrs(c[8])
            par = at.get("Parent", "")
            tgt = at.get("Target")
            if tgt:
                t = tgt.split()
                target = f"{t[0]} {t[1]} {t[2]}"
            else:
                length = (int(c[4]) - int(c[3]) + 1) // 3
                target = f"{mrna_target.get(par, par)} 1 {max(1, length)}"
            r = [c[0], a.source, "nucleotide_to_protein_match", c[3], c[4],
                 c[5] if c[5] != "." else "100", c[6], ".", f"ID={a.source}.{par};Target={target}"]
            out.write("\t".join(r) + "\n")
            n += 1
    print(f"{n} protein alignment segments")


def cmd_transcripts(a):
    tx = OrderedDict()
    for c in read_gff(a.gtf):
        if c[2] != "exon":
            continue
        at = parse_attrs(c[8])
        tid = at.get("transcript_id") or at.get("Parent")
        tx.setdefault(tid, []).append(c)
    with open(a.out, "w") as out:
        for tid, exons in tx.items():
            exons.sort(key=lambda e: int(e[3]))
            strand = exons[0][6]
            order = exons if strand != "-" else list(reversed(exons))
            pos = 1
            for e in order:
                ln = int(e[4]) - int(e[3]) + 1
                tstart, tend = pos, pos + ln - 1
                pos = tend + 1
                r = [e[0], a.source, "cDNA_match", e[3], e[4], ".", strand if strand in "+-" else "+", ".",
                     f"ID={a.source}.{tid};Target={tid} {tstart} {tend} +"]
                out.write("\t".join(r) + "\n")
    print(f"{len(tx)} transcripts")


def cmd_repeats(a):
    with open(a.out, "w") as out:
        for i, c in enumerate(read_gff(a.gff), 1):
            r = [c[0], "RepeatMasker", "repeat", c[3], c[4], ".", c[6], ".", f"ID=rep{i}"]
            out.write("\t".join(r) + "\n")


def cmd_weights(a):
    with open(a.out, "w") as out:
        for line in a.lines:
            cls, src, w = line.split(":")
            out.write(f"{cls}\t{src}\t{w}\n")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("genepred"); p.add_argument("--gff", required=True); p.add_argument("--source", required=True); p.add_argument("--out", required=True)
    p = sub.add_parser("miniprot"); p.add_argument("--gff", required=True); p.add_argument("--source", default="miniprot"); p.add_argument("--out", required=True)
    p = sub.add_parser("transcripts"); p.add_argument("--gtf", required=True); p.add_argument("--source", required=True); p.add_argument("--out", required=True)
    p = sub.add_parser("repeats"); p.add_argument("--gff", required=True); p.add_argument("--out", required=True)
    p = sub.add_parser("weights"); p.add_argument("--out", required=True); p.add_argument("lines", nargs="*", help="CLASS:SOURCE:WEIGHT")
    a = ap.parse_args()
    {"genepred": cmd_genepred, "miniprot": cmd_miniprot, "transcripts": cmd_transcripts,
     "repeats": cmd_repeats, "weights": cmd_weights}[a.cmd](a)


if __name__ == "__main__":
    main()
