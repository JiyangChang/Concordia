#!/usr/bin/env python3
"""Keep only protein-coding genes of a homolog annotation (standard library only).

  gff    clean a GFF3: keep genes of type 'gene' / 'protein_coding_gene' whose biotype (if any) is
         protein_coding and that are not flagged pseudo; keep only their transcripts with a CDS (and
         a protein_coding transcript biotype, if given). Gene / transcript types are normalised to
         'gene' / 'mRNA'; all children of kept transcripts (exon, CDS, UTR, ...) are kept unchanged.
         Writes the cleaned GFF3, the kept transcript IDs and a report (kept / dropped by type).
  fasta  restrict a user-supplied protein FASTA to the kept transcript IDs (first header word).
         If fewer than --min-match of its records match (different ID scheme), exit with code 3
         so the caller can fall back to proteins translated from the cleaned GFF3.
"""
import argparse
import os
import sys
from collections import Counter, OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import parse_attrs, read_fasta, write_fasta  # noqa: E402

GENE_OK = {"gene", "protein_coding_gene"}
BIOTYPE_KEYS = ("biotype", "gene_biotype", "gene_type", "locus_type")
TX_BIOTYPE_KEYS = ("transcript_biotype", "transcript_type", "biotype")
TX_SKIP_CHILD = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR", "start_codon", "stop_codon",
                 "UTR", "intron", "pseudogenic_exon"}


def is_gene(t):
    return t == "gene" or t.endswith("_gene") or t == "pseudogene"


def flagged_pseudo(a):
    return a.get("pseudo", "").lower() == "true" or "pseudogene" in a.get("gene_biotype", "").lower()


def biotype_ok(a, keys):
    for k in keys:
        if k in a:
            return a[k].lower() == "protein_coding"
    return True


def cmd_gff(a):
    rows = []
    with open(a.gff) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) == 9:
                rows.append((c, parse_attrs(c[8])))

    genes = OrderedDict()                 # gene ID -> (row, attrs)
    for c, at in rows:
        if is_gene(c[2]) and "ID" in at:
            genes[at["ID"]] = (c, at)
    tx = OrderedDict()                    # transcript ID -> (row, attrs, gene ID)
    for c, at in rows:
        if is_gene(c[2]) or c[2] in TX_SKIP_CHILD or "ID" not in at:
            continue
        par = at.get("Parent", "").split(",")[0]
        if par in genes:
            tx[at["ID"]] = (c, at, par)
    children = {}
    has_cds = set()
    for c, at in rows:
        for p in at.get("Parent", "").split(","):
            if p in tx and at.get("ID") not in tx:
                children.setdefault(p, []).append(c)
                if c[2] == "CDS":
                    has_cds.add(p)

    gene_tx = {}                          # gene ID -> its transcript IDs (avoids O(genes x transcripts))
    for t, (_, _, g) in tx.items():
        gene_tx.setdefault(g, []).append(t)

    if a.passthrough:                     # filter disabled: copy, list every transcript with a CDS
        with open(a.gff) as src, open(a.out_gff, "w") as dst:
            dst.write(src.read())
        with open(a.out_ids, "w") as fh:
            fh.writelines(t + "\n" for t in tx if t in has_cds)
        with open(a.out_report, "w") as fh:
            fh.write("item\tcount\nfilter disabled (homolog_protein_coding_only: false)\t0\n")
        return

    report = Counter()
    keep_genes, keep_tx = OrderedDict(), OrderedDict()
    for gid, (c, at) in genes.items():
        if c[2] not in GENE_OK:
            report[f"dropped gene: type {c[2]}"] += 1
            continue
        if not biotype_ok(at, BIOTYPE_KEYS) or flagged_pseudo(at):
            report["dropped gene: non-coding biotype / pseudo"] += 1
            continue
        kept = [t for t in gene_tx.get(gid, []) if t in has_cds
                and biotype_ok(tx[t][1], TX_BIOTYPE_KEYS) and not flagged_pseudo(tx[t][1])]
        if not kept:
            report["dropped gene: no protein-coding transcript"] += 1
            continue
        keep_genes[gid] = c
        for t in kept:
            keep_tx[t] = gid
        dropped_tx = len(gene_tx.get(gid, [])) - len(kept)
        if dropped_tx:
            report["dropped transcripts of kept genes (no CDS / non-coding)"] += dropped_tx
    report["kept genes"] = len(keep_genes)
    report["kept transcripts"] = len(keep_tx)

    with open(a.out_gff, "w") as out:
        out.write("##gff-version 3\n")
        for gid, c in keep_genes.items():
            out.write("\t".join(c[:2] + ["gene"] + c[3:]) + "\n")
            for t in (x for x in gene_tx.get(gid, []) if x in keep_tx):
                tc = tx[t][0]
                out.write("\t".join(tc[:2] + ["mRNA"] + tc[3:]) + "\n")
                for ch in children.get(t, []):
                    out.write("\t".join(ch) + "\n")
    with open(a.out_ids, "w") as fh:
        for t in keep_tx:
            fh.write(t + "\n")
    with open(a.out_report, "w") as fh:
        fh.write("item\tcount\n")
        for k in ["kept genes", "kept transcripts"] + sorted(k for k in report if k.startswith("dropped")):
            fh.write(f"{k}\t{report[k]}\n")
    print(f"kept {len(keep_genes)} genes / {len(keep_tx)} transcripts")


def cmd_fasta(a):
    ids = {l.strip() for l in open(a.ids) if l.strip()}
    recs = list(read_fasta(a.fasta))
    hit = [(h, s) for h, s in recs if h.split()[0] in ids]
    frac = len(hit) / len(recs) if recs else 0.0
    print(f"{len(hit)}/{len(recs)} protein records match kept transcripts ({100 * frac:.1f}%)")
    if frac < a.min_match:
        sys.exit(3)
    with open(a.out, "w") as fh:
        for h, s in hit:
            write_fasta(fh, h.split()[0], s.rstrip("*."))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("gff")
    p.add_argument("--gff", required=True)
    p.add_argument("--out-gff", required=True)
    p.add_argument("--out-ids", required=True)
    p.add_argument("--out-report", required=True)
    p.add_argument("--passthrough", action="store_true", help="do not filter (copy the annotation)")
    p = sub.add_parser("fasta")
    p.add_argument("--fasta", required=True)
    p.add_argument("--ids", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--min-match", type=float, default=0.5)
    a = ap.parse_args()
    {"gff": cmd_gff, "fasta": cmd_fasta}[a.cmd](a)


if __name__ == "__main__":
    main()
