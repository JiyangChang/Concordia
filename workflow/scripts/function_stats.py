#!/usr/bin/env python3
"""06_Function statistics for the report (gene level; standard library only).

Inputs: final GFF3 (gene <- mRNA mapping, gene total), InterProScan TSV chunks, eggNOG-mapper
annotations, optional GO ontology (go-basic.obo) for term names and namespaces.
Writes into --outdir:
  overview.tsv           genes annotated by InterProScan / eggNOG / both / either / neither
  interpro_databases.tsv genes per InterProScan member database
  top_pfam.tsv, top_interpro.tsv   most frequent domains / entries (with descriptions)
  go_top.tsv             most frequent GO terms (name + namespace if an OBO file is given)
  go_namespace.tsv       genes with >=1 term per namespace (only with an OBO file)
  kegg.tsv               genes with KO / pathway / EC;  kegg_top_pathways.tsv
  cog_categories.tsv     genes per COG functional category (with standard descriptions)
  eggnog_levels.tsv      genes per eggNOG max_annot_lvl
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_gff, parse_attrs  # noqa: E402

COG = {
    "J": "Translation, ribosomal structure and biogenesis", "A": "RNA processing and modification",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics", "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure", "V": "Defense mechanisms", "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility", "Z": "Cytoskeleton",
    "W": "Extracellular structures", "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons", "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism", "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only", "S": "Function unknown",
}
NAMESPACE = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}


def split(v, sep=","):
    return [x for x in (v or "").replace("|", sep).split(sep) if x and x != "-"]


def read_obo(path):
    terms, cur = {}, {}
    if not path or not os.path.exists(path):
        return terms
    for line in open(path):
        line = line.strip()
        if line == "[Term]":
            cur = {}
        elif line.startswith("id: GO:"):
            cur["id"] = line[4:]
        elif line.startswith("name: ") and "id" in cur:
            cur["name"] = line[6:]
        elif line.startswith("namespace: ") and "id" in cur:
            terms[cur["id"]] = (cur.get("name", ""), NAMESPACE.get(line[11:], line[11:]))
        elif line.startswith("alt_id: GO:") and "id" in cur:
            cur.setdefault("alt", []).append(line[8:])
            if cur["id"] in terms:
                terms[line[8:]] = terms[cur["id"]]
    return terms


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--interproscan", nargs="*", default=[])
    ap.add_argument("--eggnog", required=True)
    ap.add_argument("--go-obo", default="")
    ap.add_argument("--kegg-names", default="", help="KEGG list/pathway table: mapXXXXX<TAB>name")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    genes, tx2gene = set(), {}
    for c in read_gff(a.gff):
        at = parse_attrs(c[8])
        if c[2] == "gene":
            genes.add(at["ID"])
        elif c[2] == "mRNA":
            tx2gene[at["ID"]] = at.get("Parent", at["ID"])
    gene_of = lambda pid: tx2gene.get(pid, pid)
    n = len(genes) or 1
    pct = lambda x: f"{100 * x / n:.2f}"

    ipr_genes, member = set(), defaultdict(set)
    pfam, pfam_desc, ipr, ipr_desc = defaultdict(set), {}, defaultdict(set), {}
    go = defaultdict(set)
    for p in a.interproscan:
        if not os.path.exists(p):
            continue
        for line in open(p):
            c = line.rstrip("\n").split("\t")
            if len(c) < 11:
                continue
            g = gene_of(c[0])
            ipr_genes.add(g)
            member[c[3]].add(g)
            if c[3] == "Pfam":
                pfam[c[4]].add(g)
                pfam_desc.setdefault(c[4], c[5])
            if len(c) > 11 and c[11] not in ("", "-"):
                ipr[c[11]].add(g)
                ipr_desc.setdefault(c[11], c[12] if len(c) > 12 else "")
            if len(c) > 13:
                for t in split(c[13]):
                    go[t.split("(")[0]].add(g)

    egg_genes, ko, pathway, ec, cog, lvl = set(), set(), defaultdict(set), set(), defaultdict(set), Counter()
    header = None
    if os.path.exists(a.eggnog):
        for line in open(a.eggnog):
            if line.startswith("#query"):
                header = line[1:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or header is None:
                continue
            d = dict(zip(header, line.rstrip("\n").split("\t")))
            g = gene_of(d["query"])
            egg_genes.add(g)
            for t in split(d.get("GOs")):
                go[t].add(g)
            if split(d.get("KEGG_ko")):
                ko.add(g)
            for pw in split(d.get("KEGG_Pathway")):
                if pw.startswith("map"):
                    pathway[pw].add(g)
            if split(d.get("EC")):
                ec.add(g)
            for ch in (d.get("COG_category") or "").replace("-", ""):
                cog[ch].add(g)
            if d.get("max_annot_lvl") not in (None, "", "-"):
                lvl[d["max_annot_lvl"]] += 1

    both = ipr_genes & egg_genes
    either = ipr_genes | egg_genes
    write(os.path.join(a.outdir, "overview.tsv"), ["category", "genes", "percent"], [
        ["total genes", len(genes), "100.00"],
        ["InterProScan (any member database)", len(ipr_genes), pct(len(ipr_genes))],
        ["eggNOG-mapper", len(egg_genes), pct(len(egg_genes))],
        ["both", len(both), pct(len(both))],
        ["InterProScan only", len(ipr_genes - egg_genes), pct(len(ipr_genes - egg_genes))],
        ["eggNOG-mapper only", len(egg_genes - ipr_genes), pct(len(egg_genes - ipr_genes))],
        ["either (annotated)", len(either), pct(len(either))],
        ["neither (unannotated)", len(genes - either), pct(len(genes - either))],
        ["with GO term(s)", len(set().union(*go.values())) if go else 0,
         pct(len(set().union(*go.values()))) if go else "0.00"],
    ])
    write(os.path.join(a.outdir, "interpro_databases.tsv"), ["database", "genes", "percent"],
          [[db, len(gs), pct(len(gs))] for db, gs in sorted(member.items(), key=lambda x: -len(x[1]))])
    write(os.path.join(a.outdir, "top_pfam.tsv"), ["pfam", "description", "genes"],
          [[k, pfam_desc.get(k, ""), len(v)] for k, v in sorted(pfam.items(), key=lambda x: -len(x[1]))[:a.top]])
    write(os.path.join(a.outdir, "top_interpro.tsv"), ["interpro", "description", "genes"],
          [[k, ipr_desc.get(k, ""), len(v)] for k, v in sorted(ipr.items(), key=lambda x: -len(x[1]))[:a.top]])

    obo = read_obo(a.go_obo)
    write(os.path.join(a.outdir, "go_top.tsv"), ["go", "name", "namespace", "genes"],
          [[t, obo.get(t, ("-", "-"))[0], obo.get(t, ("-", "-"))[1], len(gs)]
           for t, gs in sorted(go.items(), key=lambda x: -len(x[1]))[:a.top]])
    ns_rows = []
    if obo:
        by_ns = defaultdict(set)
        for t, gs in go.items():
            if t in obo:
                by_ns[obo[t][1]] |= gs
        ns_rows = [[k, len(by_ns[k]), pct(len(by_ns[k]))] for k in ("BP", "MF", "CC")]
    write(os.path.join(a.outdir, "go_namespace.tsv"), ["namespace", "genes", "percent"], ns_rows)

    pw_genes = set().union(*pathway.values()) if pathway else set()
    write(os.path.join(a.outdir, "kegg.tsv"), ["category", "genes", "percent"], [
        ["with KEGG KO", len(ko), pct(len(ko))],
        ["with KEGG pathway", len(pw_genes), pct(len(pw_genes))],
        ["with EC number", len(ec), pct(len(ec))]])
    kegg_names = {}
    if a.kegg_names and os.path.exists(a.kegg_names):
        for line in open(a.kegg_names):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                kegg_names[p[0].replace("path:", "")] = p[1]
    write(os.path.join(a.outdir, "kegg_top_pathways.tsv"), ["pathway", "name", "genes"],
          [[k, kegg_names.get(k, "-"), len(v)] for k, v in sorted(pathway.items(), key=lambda x: -len(x[1]))[:a.top]])
    write(os.path.join(a.outdir, "cog_categories.tsv"), ["code", "description", "genes", "percent"],
          [[k, COG.get(k, "-"), len(v), pct(len(v))] for k, v in sorted(cog.items(), key=lambda x: -len(x[1]))])
    write(os.path.join(a.outdir, "eggnog_levels.tsv"), ["max_annot_lvl", "genes"], lvl.most_common(a.top))


if __name__ == "__main__":
    main()
