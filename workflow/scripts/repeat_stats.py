#!/usr/bin/env python3
"""01_Repeat statistics for the report (standard library only).

Inputs (all optional except --te-gff and --tbl):
  --te-gff     EDTA whole-genome TE annotation (GFF3 with Classification= / Method= attributes)
  --tbl        RepeatMasker .tbl (genome size, masked bp, simple repeats, low complexity)
  --lib        EDTA TE library FASTA (headers name#Class/Superfamily)
  --divsum     RepeatMasker calcDivergenceFromAlign .divsum (Kimura repeat landscape)
  --pass-list  LTR_retriever / EDTA intact-LTR pass.list (insertion times)
Writes into --outdir:
  overview.tsv    genome size, masked bp/% (RepeatMasker), TE bp/% (EDTA), simple / low complexity
  te_classes.tsv  per class and superfamily: copies, merged bp, % of genome
  te_methods.tsv  EDTA structural (intact) vs homology annotations: copies, merged bp
  library.tsv     consensus sequences per class in the TE library
  landscape.tsv   divergence bin x class -> bp and % of genome
  ltr_age.tsv     intact LTRs: superfamily, identity, insertion time (Mya)
"""
import argparse
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_gff, parse_attrs  # noqa: E402

SKIP_TYPES = {"repeat_region", "long_terminal_repeat", "target_site_duplication"}


def merged_bp(ivs):
    """Total bp covered by (seqid, start, end) intervals, overlaps counted once."""
    total, cur = 0, None
    for seq, s, e in sorted(ivs):
        if cur and seq == cur[0] and s <= cur[2] + 1:
            cur[2] = max(cur[2], e)
        else:
            if cur:
                total += cur[2] - cur[1] + 1
            cur = [seq, s, e]
    if cur:
        total += cur[2] - cur[1] + 1
    return total


def parse_tbl(path):
    d = {}
    if not path or not os.path.exists(path):
        return d
    for line in open(path):
        m = re.match(r"total length:\s+(\d+) bp", line.strip())
        if m:
            d["genome_bp"] = int(m.group(1))
        m = re.match(r"bases masked:\s+(\d+) bp\s+\(\s*([\d.]+) %\)", line.strip())
        if m:
            d["masked_bp"], d["masked_pct"] = int(m.group(1)), m.group(2)
        m = re.match(r"(Simple repeats|Low complexity):\s+(\d+)\s+(\d+) bp\s+([\d.]+) %", line.strip())
        if m:
            d[m.group(1)] = (m.group(2), m.group(3), m.group(4))
    return d


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--te-gff", required=True)
    ap.add_argument("--tbl", required=True)
    ap.add_argument("--lib", default="")
    ap.add_argument("--divsum", default="")
    ap.add_argument("--pass-list", default="")
    ap.add_argument("--mutation-rate", default="1.3e-8 (EDTA default)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    tbl = parse_tbl(a.tbl)
    gsize = tbl.get("genome_bp", 0)
    pct = lambda bp: f"{100 * bp / gsize:.2f}" if gsize else "-"

    # ---- EDTA TE annotation
    by_cls, by_sf, by_method = defaultdict(list), defaultdict(list), defaultdict(list)
    n_cls, n_sf, n_method = Counter(), Counter(), Counter()
    all_iv = []
    for c in read_gff(a.te_gff):
        if c[2] in SKIP_TYPES:
            continue
        at = parse_attrs(c[8])
        cl = at.get("Classification") or at.get("classification") or c[2]
        top = cl.split("/")[0]
        iv = (c[0], int(c[3]), int(c[4]))
        by_cls[top].append(iv)
        by_sf[cl].append(iv)
        n_cls[top] += 1
        n_sf[cl] += 1
        m = at.get("Method") or at.get("method") or "unknown"
        by_method[m].append(iv)
        n_method[m] += 1
        all_iv.append(iv)
    te_bp = merged_bp(all_iv)

    rows = []
    for top in sorted(by_cls, key=lambda k: -merged_bp(by_cls[k])):
        bp = merged_bp(by_cls[top])
        rows.append([top, "(all)", n_cls[top], bp, pct(bp)])
        for sf in sorted((k for k in by_sf if k.split("/")[0] == top), key=lambda k: -merged_bp(by_sf[k])):
            if sf == top:
                continue
            bp = merged_bp(by_sf[sf])
            rows.append([top, sf.split("/", 1)[1] if "/" in sf else sf, n_sf[sf], bp, pct(bp)])
    rows.append(["total TE", "-", sum(n_cls.values()), te_bp, pct(te_bp)])
    write(os.path.join(a.outdir, "te_classes.tsv"), ["class", "superfamily", "copies", "bp", "genome_%"], rows)
    write(os.path.join(a.outdir, "te_methods.tsv"), ["method", "copies", "bp", "genome_%"],
          [[m, n_method[m], merged_bp(v), pct(merged_bp(v))] for m, v in sorted(by_method.items())])

    ov = [["genome size (bp)", gsize or "-", ""],
          ["masked by RepeatMasker (EDTA library, incl. simple/low complexity)", tbl.get("masked_bp", "-"),
           tbl.get("masked_pct", "-")],
          ["TE annotated by EDTA (merged)", te_bp, pct(te_bp)]]
    for k in ("Simple repeats", "Low complexity"):
        if k in tbl:
            ov.append([f"{k} ({tbl[k][0]} elements)", tbl[k][1], tbl[k][2]])
    write(os.path.join(a.outdir, "overview.tsv"), ["item", "bp", "genome_%"], ov)

    # ---- TE library
    lib = Counter()
    if a.lib and os.path.exists(a.lib):
        for line in open(a.lib):
            if line.startswith(">"):
                cl = line[1:].split()[0].split("#", 1)[1] if "#" in line else "unknown"
                lib[cl] += 1
    write(os.path.join(a.outdir, "library.tsv"), ["classification", "consensus_sequences"],
          sorted(lib.items(), key=lambda x: -x[1]))

    # ---- Kimura landscape
    land = []
    if a.divsum and os.path.exists(a.divsum):
        lines = open(a.divsum).read().splitlines()
        try:
            i = next(j for j, l in enumerate(lines) if l.startswith("Div "))
            header = lines[i].split()[1:]
            for l in lines[i + 1:]:
                p = l.split()
                if not p or not p[0].isdigit():
                    break
                agg = Counter()
                for name, v in zip(header, p[1:]):
                    agg[name.split("/")[0]] += int(float(v))
                for cl, bp in sorted(agg.items()):
                    land.append([p[0], cl, bp, pct(bp)])
        except StopIteration:
            pass
    write(os.path.join(a.outdir, "landscape.tsv"), ["kimura_div", "class", "bp", "genome_%"], land)

    # ---- intact LTR insertion times
    ages = []
    if a.pass_list and os.path.exists(a.pass_list):
        header = None
        for line in open(a.pass_list):
            p = line.split()                      # LTR_retriever writes the header space-separated
            if line.startswith("#"):
                header = [h.lstrip("#").strip() for h in p]
                continue
            if header is None or len(p) < len(header):
                continue
            d = dict(zip(header, p))
            t = d.get("Insertion_Time") or d.get("Insertion_time")
            try:
                mya = float(t) / 1e6
            except (TypeError, ValueError):
                continue
            ages.append([d.get("LTR_loc", ""), d.get("SuperFamily", "unknown"), d.get("Identity", ""), f"{mya:.3f}"])
    with open(os.path.join(a.outdir, "ltr_age.tsv"), "w") as fh:
        fh.write(f"# mutation rate u = {a.mutation_rate}; T = K / (2u)\n")
        fh.write("ltr\tsuperfamily\tidentity\tinsertion_mya\n")
        for r in ages:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
