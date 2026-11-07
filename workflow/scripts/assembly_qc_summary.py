#!/usr/bin/env python3
"""Collect all genome assembly quality metrics into one table: section, metric, value.

Every input is optional; missing ones are reported as 'not run'.
"""
import argparse
import glob
import os
import re

BUSCO_RE = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)")


def rows_kv(path):
    if path and os.path.exists(path):
        return [l.rstrip("\n").split("\t") for l in open(path) if "\t" in l]
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True)
    ap.add_argument("--busco", nargs="*", default=[], help="lineage=short_summary.txt")
    ap.add_argument("--lai", default="")
    ap.add_argument("--merqury-qv", default="")
    ap.add_argument("--merqury-completeness", default="")
    ap.add_argument("--telomere", default="")
    ap.add_argument("--fcs", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = []

    stats = rows_kv(a.stats)
    for r in stats[1:]:
        if r[0] in ("total_length", "count", "N50", "L50", "N90", "GC_percent", "N_bases"):
            out.append(("contiguity", f"scaffold_{r[0]}", r[1]))
            out.append(("contiguity", f"contig_{r[0]}", r[2]))
        elif r[0] in ("gaps", "gap_length"):
            out.append(("contiguity", r[0], r[1]))

    for item in a.busco:
        lin, path = item.split("=", 1)
        m = BUSCO_RE.search(open(path).read()) if os.path.exists(path) else None
        out.append(("completeness", f"BUSCO_genome_{lin}",
                    f"C:{m.group(1)}% [S:{m.group(2)}%, D:{m.group(3)}%], F:{m.group(4)}%, M:{m.group(5)}%, n:{m.group(6)}"
                    if m else "not run"))

    if a.lai:  # LAI runs after EDTA; "not applicable" when LTR content is too low
        lai = "no whole_genome line"
        if os.path.exists(a.lai):
            for line in open(a.lai):
                c = line.split()
                if line.startswith("# LAI not applicable"):
                    m = re.search(r"content \(([\d.]+)%\)", line)
                    lai = (f"not applicable (LTR content {m.group(1)}% < 5%)" if m else
                           "not applicable: " + line.split(":", 1)[1].strip()[:80])
                elif c and c[0] == "whole_genome":
                    lai = f"raw_LAI={c[-2]}, LAI={c[-1]}"
        out.append(("completeness", "LAI", lai))

    qv = "not run"
    if a.merqury_qv and os.path.exists(a.merqury_qv):
        c = open(a.merqury_qv).readline().split()
        if len(c) >= 4:
            qv = f"QV={c[3]} (error rate {c[4] if len(c) > 4 else '-'})"
    out.append(("accuracy", "Merqury_QV", qv))
    comp = "not run"
    if a.merqury_completeness and os.path.exists(a.merqury_completeness):
        c = open(a.merqury_completeness).readline().split()
        if len(c) >= 5:
            comp = f"{c[4]}%"
    out.append(("completeness", "Merqury_kmer_completeness", comp))

    tel = dict(rows_kv(a.telomere)) if a.telomere else {}
    if tel:
        out.append(("telomere", "motif", tel.get("motif", "")))
        out.append(("telomere", "sequences_with_telomeres_at_both_ends", tel.get("telomere_both_ends", "")))
        out.append(("telomere", "sequences_with_telomere_at_one_end", tel.get("telomere_one_end", "")))
    else:
        out.append(("telomere", "tidk", "not run"))

    if a.fcs and os.path.exists(a.fcs):
        actions = {}
        for line in open(a.fcs):
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) > 4:
                actions[c[4]] = actions.get(c[4], 0) + 1
        out.append(("contamination", "FCS-GX_flagged_regions",
                    ", ".join(f"{k}={v}" for k, v in sorted(actions.items())) or "none"))
    else:
        out.append(("contamination", "FCS-GX", "not run"))

    with open(a.out, "w") as fh:
        fh.write("section\tmetric\tvalue\n")
        for r in out:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
