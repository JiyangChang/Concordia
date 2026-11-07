#!/usr/bin/env python3
"""06_Function: merge InterProScan + eggNOG-mapper results and write them back to the GFF3."""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import parse_attrs  # noqa: E402

FIELDS = ["InterPro", "Pfam", "GO", "KEGG_ko", "KEGG_Pathway", "EC", "COG_category", "Preferred_name", "Description"]
NAME_FIELDS = {"Preferred_name", "Description"}   # gene_name= / product= in the GFF3, not annotation categories


def _split(v, sep=","):
    return [x for x in (v or "").replace("|", sep).split(sep) if x and x != "-"]


def load_interproscan(paths, ann):
    for p in paths:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            continue
        with open(p) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 11:
                    continue
                pid = c[0]
                if c[3] == "Pfam":
                    ann[pid]["Pfam"].add(c[4])
                if len(c) > 11 and c[11] not in ("", "-"):
                    ann[pid]["InterPro"].add(c[11])
                if len(c) > 13:
                    for go in _split(c[13]):
                        ann[pid]["GO"].add(go.split("(")[0])


def load_eggnog(path, ann):
    if not os.path.exists(path):
        return
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#query"):
                header = line[1:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or header is None:
                continue
            d = dict(zip(header, line.rstrip("\n").split("\t")))
            pid = d["query"]
            for go in _split(d.get("GOs")):
                ann[pid]["GO"].add(go)
            for ko in _split(d.get("KEGG_ko")):
                ann[pid]["KEGG_ko"].add(ko.replace("ko:", ""))
            for pw in _split(d.get("KEGG_Pathway")):
                if pw.startswith("map"):
                    ann[pid]["KEGG_Pathway"].add(pw)
            for ec in _split(d.get("EC")):
                ann[pid]["EC"].add(ec)
            for cog in (d.get("COG_category") or "").replace("-", ""):
                ann[pid]["COG_category"].add(cog)
            if d.get("Preferred_name") not in (None, "", "-"):
                ann[pid]["Preferred_name"].add(d["Preferred_name"])
            if d.get("Description") not in (None, "", "-"):
                ann[pid]["Description"].add(d["Description"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--eggnog", required=True)
    ap.add_argument("--interproscan", nargs="*", default=[])
    ap.add_argument("--out-gff", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-go", required=True)
    ap.add_argument("--out-kegg", required=True)
    ap.add_argument("--out-stats", required=True)
    a = ap.parse_args()

    ann = defaultdict(lambda: defaultdict(set))
    load_interproscan(a.interproscan, ann)
    load_eggnog(a.eggnog, ann)

    mrna2gene, genes = {}, []
    with open(a.gff) as fh, open(a.out_gff, "w") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                out.write(line); continue
            c = line.rstrip("\n").split("\t")
            if c[2] == "gene":
                genes.append(parse_attrs(c[8])["ID"])
            if c[2] == "mRNA":
                at = parse_attrs(c[8])
                mrna2gene[at["ID"]] = at.get("Parent", "")
                d = ann.get(at["ID"])
                if d:
                    dbx = [f"InterPro:{x}" for x in sorted(d["InterPro"])] + [f"Pfam:{x}" for x in sorted(d["Pfam"])]
                    if dbx:
                        c[8] += ";Dbxref=" + ",".join(dbx)
                    if d["GO"]:
                        c[8] += ";Ontology_term=" + ",".join(sorted(d["GO"]))
                    if d["Preferred_name"]:
                        c[8] += ";gene_name=" + sorted(d["Preferred_name"])[0].replace(";", ",")
                    if d["Description"]:
                        c[8] += ";product=" + sorted(d["Description"])[0].replace(";", ",").replace("=", " ")
            out.write("\t".join(c) + "\n")

    with open(a.out_tsv, "w") as t, open(a.out_go, "w") as g, open(a.out_kegg, "w") as k:
        t.write("transcript\tgene\t" + "\t".join(FIELDS) + "\n")
        for tid, gid in mrna2gene.items():
            d = ann.get(tid, {})
            t.write(f"{tid}\t{gid}\t" + "\t".join(",".join(sorted(d.get(f, []))) or "-" for f in FIELDS) + "\n")
            for go in sorted(d.get("GO", [])):
                g.write(f"{gid}\t{go}\n")
            for ko in sorted(d.get("KEGG_ko", [])):
                k.write(f"{gid}\t{ko}\n")

    gene_any = defaultdict(set)
    for tid, gid in mrna2gene.items():
        for f in FIELDS:
            if ann.get(tid, {}).get(f):
                gene_any[f].add(gid)
                if f not in NAME_FIELDS:          # "any" = any annotation category, not just a name/description
                    gene_any["any"].add(gid)
    n = len(genes) or 1
    with open(a.out_stats, "w") as s:
        s.write("database\tgenes\tpercent\n")
        for f in ["any"] + FIELDS:
            s.write(f"{f}\t{len(gene_any[f])}\t{100 * len(gene_any[f]) / n:.2f}\n")


if __name__ == "__main__":
    main()
