#!/usr/bin/env python3
"""05_EVM/final: remove TE-derived genes and too-short models from the refined gene set.

A gene is removed for TE reasons only when BOTH hold:
  1. TE overlap > --max-te-overlap: fraction of the gene's coding bases (union of CDS of all its
     transcripts) covered by EDTA TE annotation. CDS bases are used instead of the gene span, so
     TEs that sit only in introns or UTRs do not count.
  2. no RNA evidence. A gene has RNA evidence if ANY of these is true:
     a. StringTie: one of its exons overlaps an exon of a StringTie transcript (same strand, or
        the transcript is unstranded)
     b. EviAnn: refine marked the gene eviann_support=yes (added from EviAnn, merged with an EviAnn
        model of the same CDS intron chain, or given EviAnn isoforms); older refine tables with
        gffcompare '=' matches are still accepted
     c. BRAKER4: its CDS overlaps (same strand) a BRAKER4 model that gene_support.tsv reports as
        RNA-supported
With --protein-evidence (annotation without RNA-seq), criterion 2 becomes "no RNA and no protein
evidence": protein evidence = a CDS of the gene overlaps (same strand) a CDS of a miniprot alignment or
a GeMoMa model. Without it every TE-overlapping gene would be removed, since there is no RNA evidence.
Independently, genes whose longest CDS is shorter than --min-cds-aa are removed - unless --keep-short-supported
and the gene has RNA evidence (as above; without RNA-seq: protein evidence): such genes are kept and marked
short_cds=yes (a server run removed 1,547 RNA-supported genes with min_cds_aa 100).
TEsorter protein hits are recorded in the decision table only; they never cause removal.

Outputs: filtered GFF3, removed_genes.tsv (gene, reason) and filter_decisions.tsv (every gene).
"""
import argparse
import bisect
import csv
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_gff, parse_attrs  # noqa: E402
from gff_models import load_models  # noqa: E402
from refine import write_gene, cds_len  # noqa: E402

TRUE_WORDS = {"yes", "true", "y", "t", "supported", "full", "partial"}
RNA_COL = re.compile(r"rna|isoseq|iso_seq|transcript|intron|est", re.I)


# ---------- intervals ----------
def merge(ivs):
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


class IntervalIndex:
    """Merged, sorted intervals per (seqid[, strand]) for fast coverage / overlap queries."""

    def __init__(self, items):
        groups = defaultdict(list)
        for key, s, e in items:
            groups[key].append((s, e))
        self.iv = {k: merge(v) for k, v in groups.items()}
        self.starts = {k: [x[0] for x in v] for k, v in self.iv.items()}

    def covered(self, key, s, e):
        """Number of bases of [s, e] covered by the index."""
        iv, st = self.iv.get(key), self.starts.get(key)
        if not iv:
            return 0
        i = max(0, bisect.bisect_right(st, e) - 1)
        n = 0
        while i >= 0 and iv[i][1] >= s:
            n += max(0, min(e, iv[i][1]) - max(s, iv[i][0]) + 1)
            i -= 1
        return n

    def overlaps(self, key, s, e):
        return self.covered(key, s, e) > 0


# ---------- evidence loaders ----------
def te_index(path):
    items = [(c[0], int(c[3]), int(c[4])) for c in read_gff(path)] if path and os.path.exists(path) else []
    return IntervalIndex(items)


def stringtie_index(path):
    items = []
    if path and os.path.exists(path):
        for c in read_gff(path):
            if c[2] == "exon":
                items.append(((c[0], c[6] if c[6] in "+-" else "."), int(c[3]), int(c[4])))
    return IntervalIndex(items)


def protein_index(paths):
    items = []
    for p in paths or []:
        if p and os.path.exists(p):
            items += [((c[0], c[6]), int(c[3]), int(c[4])) for c in read_gff(p) if c[2] == "CDS"]
    return IntervalIndex(items)


def eviann_matched(decisions):
    """EVM transcript and gene IDs with an exact ('=') EviAnn match."""
    ids = set()
    if decisions and os.path.exists(decisions):
        with open(decisions) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("class_code") == "=":
                    ids.update(x for x in (r.get("evm_gene"), r.get("evm_tx")) if x and x != "-")
    return ids


def _truthy(v):
    v = (v or "").strip().lower()
    if v in TRUE_WORDS:
        return True
    try:
        return float(v.rstrip("%")) > 0
    except ValueError:
        return False


def read_gene_support(path):
    """(header, rows) of BRAKER4 gene_support.tsv: '#' summary lines are skipped; the header is the
    first remaining line (a leading '#' on it is tolerated if it names gene_id)."""
    header, rows = None, []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                if header is None and "gene_id" in line:
                    header = [h.strip() for h in line.lstrip("#").rstrip("\n").split("\t")]
                continue
            p = line.rstrip("\n").split("\t")
            if header is None:
                header = [h.strip() for h in p]
            else:
                rows.append(p)
    return header or [], rows


def braker4_supported(support_path, gff_path):
    """IntervalIndex of CDS of BRAKER4 models reported as RNA-supported, plus a status note.

    BRAKER4 gene_support.tsv (scripts/gene_support_summary.py): one row per transcript with
    gene_id, transcript_id, num_introns, introns_sup_rnaseq, exons_sup_rnaseq, ... A model counts as
    RNA-supported when >= 1 intron is supported by RNA-seq hints, or (single-exon models) >= 1 exon
    is. For other layouts, every column whose name matches RNA_COL is used (value > 0 / yes).
    """
    if not support_path or not os.path.exists(support_path) or os.path.getsize(support_path) == 0:
        return IntervalIndex([]), "gene_support.tsv missing or empty: criterion c not used"
    header, rows = read_gene_support(support_path)
    ix = {h: i for i, h in enumerate(header)}
    get = lambda r, k: r[ix[k]] if k in ix and ix[k] < len(r) else ""
    ok = set()
    if {"introns_sup_rnaseq", "exons_sup_rnaseq", "num_introns"} <= set(ix):
        for r in rows:
            n_intr = get(r, "num_introns")
            sup = _truthy(get(r, "introns_sup_rnaseq")) if n_intr not in ("", "0") else _truthy(get(r, "exons_sup_rnaseq"))
            if sup:
                ok.update(x for x in (get(r, "gene_id"), get(r, "transcript_id")) if x)
        used = "introns_sup_rnaseq (single-exon: exons_sup_rnaseq)"
    else:
        cols = [i for i, h in enumerate(header) if i > 0 and RNA_COL.search(h)]
        if not cols:
            return IntervalIndex([]), f"no RNA column found in gene_support.tsv header {header}: criterion c not used"
        ok = {r[0] for r in rows if any(i < len(r) and _truthy(r[i]) for i in cols)}
        used = ", ".join(header[i] for i in cols)
    genes, _ = load_models(gff_path) if gff_path and os.path.exists(gff_path) else ({}, [])
    items = []
    for gid, g in genes.items():
        for tid, t in g["tx"].items():
            if gid in ok or tid in ok:
                items += [((f[0], f[6]), int(f[3]), int(f[4])) for f in t["feats"] if f[2] == "CDS"]
    return IntervalIndex(items), f"criterion c uses gene_support.tsv column(s): {used}; {len(ok)} supported IDs"


def tesorter_hits(path):
    hits = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                hits[c[0].split("|")[0]] = "/".join(x for x in c[1:3] if x)
    return hits


# ---------- per-gene evaluation ----------
def evaluate(gid, g, te, st, evm_eq, b4):
    cds = [(f[0], f[6], int(f[3]), int(f[4])) for t in g["tx"].values() for f in t["feats"] if f[2] == "CDS"]
    exons = [(f[0], f[6], int(f[3]), int(f[4])) for t in g["tx"].values() for f in t["feats"] if f[2] == "exon"]
    exons = exons or cds
    seqid, strand = g["row"][0], g["row"][6]
    cds_union = merge([(s, e) for _, _, s, e in cds])
    cds_bp = sum(e - s + 1 for s, e in cds_union)
    te_bp = sum(te.covered(seqid, s, e) for s, e in cds_union)

    rna_st = any(st.overlaps((c, sd), s, e) or st.overlaps((c, "."), s, e) for c, sd, s, e in exons)
    rna_evi = (g["row"][1] == "EviAnn_refine" or parse_attrs(g["row"][8]).get("eviann_support") == "yes"
               or gid in evm_eq or any(t in evm_eq for t in g["tx"]))
    rna_b4 = any(b4.overlaps((seqid, strand), s, e) for s, e in cds_union)
    return {
        "cds_bp": cds_bp,
        "longest_cds_aa": max((cds_len(t) for t in g["tx"].values()), default=0) // 3,
        "te_bp": te_bp,
        "te_frac": te_bp / cds_bp if cds_bp else 0.0,
        "rna_stringtie": rna_st,
        "rna_eviann": rna_evi,
        "rna_braker4": rna_b4,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="refined GFF3")
    ap.add_argument("--te-gff", required=True, help="EDTA TE annotation GFF3")
    ap.add_argument("--max-te-overlap", type=float, default=0.5)
    ap.add_argument("--transcripts", default="", help="StringTie GTF")
    ap.add_argument("--refine-decisions", default="", help="refine_decisions.tsv")
    ap.add_argument("--braker4-gff", default="")
    ap.add_argument("--braker4-support", default="", help="BRAKER4 gene_support.tsv (may be empty)")
    ap.add_argument("--tesorter", default="", help="TEsorter *.cls.tsv (recorded only)")
    ap.add_argument("--min-cds-aa", type=int, default=50)
    ap.add_argument("--keep-short-supported", action="store_true",
                    help="keep genes with a CDS below --min-cds-aa when they have RNA (or, without RNA-seq, protein) "
                         "evidence; marked short_cds=yes")
    ap.add_argument("--protein-evidence", nargs="*", default=None,
                    help="CDS GFF3s (miniprot, GeMoMa): protein evidence also protects TE-overlapping genes "
                         "(used when there is no RNA-seq)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--removed", required=True)
    ap.add_argument("--decisions", required=True)
    a = ap.parse_args()

    te = te_index(a.te_gff)
    st = stringtie_index(a.transcripts)
    evm_eq = eviann_matched(a.refine_decisions)
    b4, b4_note = braker4_supported(a.braker4_support, a.braker4_gff)
    tes = tesorter_hits(a.tesorter)
    use_prot = a.protein_evidence is not None
    prot = protein_index(a.protein_evidence) if use_prot else None
    print(b4_note)

    genes, _ = load_models(a.gff)
    kept = 0
    removed = defaultdict(int)
    cols = ["gene", "source", "cds_bp", "longest_cds_aa", "te_overlap_bp", "te_overlap_frac", "tesorter",
            "rna_stringtie", "rna_eviann", "rna_braker4", "rna_supported"] + (["protein"] if use_prot else []) + [
            "action", "reason"]
    with open(a.out, "w") as out, open(a.removed, "w") as rm, open(a.decisions, "w") as dec:
        out.write("##gff-version 3\n")
        rm.write("gene\treason\n")
        dec.write(f"# max_te_overlap={a.max_te_overlap} min_cds_aa={a.min_cds_aa} "
                  f"keep_short_supported={'yes' if a.keep_short_supported else 'no'}; {b4_note}"
                  + ("; no RNA-seq: protein evidence (miniprot / GeMoMa) also protects TE-overlapping genes"
                     if use_prot else "") + "\n")
        dec.write("\t".join(cols) + "\n")
        for gid, g in genes.items():
            v = evaluate(gid, g, te, st, evm_eq, b4)
            rna = v["rna_stringtie"] or v["rna_eviann"] or v["rna_braker4"]
            has_prot = use_prot and any(
                prot.overlaps((f[0], f[6]), int(f[3]), int(f[4]))
                for t in g["tx"].values() for f in t["feats"] if f[2] == "CDS")
            hit = next((tes[t] for t in g["tx"] if t in tes), "")
            if v["te_frac"] > a.max_te_overlap and not rna and not has_prot:
                action, reason = "remove", (f"TE_overlap_{v['te_frac']:.2f}_no_RNA" + ("_no_protein" if use_prot else ""))
            elif v["longest_cds_aa"] < a.min_cds_aa and a.keep_short_supported and (rna or has_prot):
                action, reason = "keep", f"short_CDS_{v['longest_cds_aa']}aa_supported"
                g["row"][8] = g["row"][8].rstrip(";") + ";short_cds=yes"
            elif v["longest_cds_aa"] < a.min_cds_aa:
                action, reason = "remove", f"CDS_{v['longest_cds_aa']}aa_lt_{a.min_cds_aa}"
            else:
                action, reason = "keep", ""
            yn = lambda x: "yes" if x else "no"
            dec.write("\t".join(map(str, [
                gid, g["row"][1], v["cds_bp"], v["longest_cds_aa"], v["te_bp"], f"{v['te_frac']:.3f}",
                hit or "-", yn(v["rna_stringtie"]), yn(v["rna_eviann"]), yn(v["rna_braker4"]), yn(rna)] +
                ([yn(has_prot)] if use_prot else []) + [action, reason or "-"])) + "\n")
            if action == "remove":
                rm.write(f"{gid}\t{reason}\n")
                removed[reason.split("_")[0]] += 1
            else:
                write_gene(out, g)
                kept += 1
    print(f"kept {kept} genes; removed: " + (", ".join(f"{k}={n}" for k, n in removed.items()) or "none"))


if __name__ == "__main__":
    main()
