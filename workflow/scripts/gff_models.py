"""Group a GFF3 into gene -> transcript -> features (standard library only)."""
import os
import sys
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_gff, parse_attrs  # noqa: E402

TX_TYPES = {"mRNA", "transcript", "lnc_RNA", "lncRNA", "ncRNA", "pseudogenic_transcript"}


def is_gene_type(t):
    """'gene' plus the typed variants used by e.g. Araport11 / IRGSP (protein_coding_gene, ncRNA_gene,
    transposable_element_gene, pseudogene)."""
    return t == "gene" or t.endswith("_gene") or t == "pseudogene"


def load_models(path):
    """Return OrderedDict gene_id -> {"row": cols, "tx": OrderedDict tx_id -> {"row", "feats": []}}."""
    genes, tx_parent, orphans = OrderedDict(), {}, []
    rows = list(read_gff(path))
    for c in rows:
        a = parse_attrs(c[8])
        if is_gene_type(c[2]):
            genes[a["ID"]] = {"row": c, "tx": OrderedDict()}
    for c in rows:
        a = parse_attrs(c[8])
        if c[2] in TX_TYPES:
            par = a.get("Parent", "").split(",")[0]
            tid = a["ID"]
            if par not in genes:  # transcript without gene: create one
                genes[par or tid + ".gene"] = {"row": [c[0], c[1], "gene", c[3], c[4], c[5], c[6], ".",
                                                     f"ID={par or tid + '.gene'}"], "tx": OrderedDict()}
                par = par or tid + ".gene"
            genes[par]["tx"][tid] = {"row": c, "feats": []}
            tx_parent[tid] = par
    for c in rows:
        if is_gene_type(c[2]) or c[2] in TX_TYPES:
            continue
        a = parse_attrs(c[8])
        for p in a.get("Parent", "").split(","):
            if p in tx_parent:
                genes[tx_parent[p]]["tx"][p]["feats"].append(c)
            else:
                orphans.append(c)
    return genes, orphans


def is_coding(gene):
    return any(f[2] == "CDS" for t in gene["tx"].values() for f in t["feats"])
