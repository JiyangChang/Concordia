#!/usr/bin/env python3
"""Summarise tidk search windows: does each sequence have telomeric repeats at its start / end?

A sequence end is called telomeric when the first (or last) window holds >= --min-repeats motif
copies on either strand. This is a heuristic; the per-sequence table lets the user check it.
"""
import argparse
import csv
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, help="tidk search *_telomeric_repeat_windows.tsv")
    ap.add_argument("--seqs", required=True, help="assembly_stats per-sequence table (seq_id, length)")
    ap.add_argument("--min-repeats", type=int, default=50)
    ap.add_argument("--motif", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-summary", required=True)
    a = ap.parse_args()

    lengths = {}
    with open(a.seqs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lengths[r["seq_id"]] = int(r["length"])

    win = defaultdict(list)
    with open(a.windows) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = [h.strip().lower() for h in rows[0]]
    i_id = 0
    i_win = next(i for i, h in enumerate(header) if "window" in h)
    i_cnt = [i for i, h in enumerate(header) if "forward" in h or "reverse" in h]
    for r in rows[1:]:
        if len(r) <= max([i_win] + i_cnt):
            continue
        win[r[i_id]].append((int(float(r[i_win])), max(int(float(r[i])) for i in i_cnt)))

    n_both = n_one = 0
    with open(a.out, "w") as fh:
        fh.write("seq_id\tlength\tstart_repeats\tend_repeats\ttelomere_start\ttelomere_end\n")
        for sid, length in sorted(lengths.items(), key=lambda x: -x[1]):
            w = sorted(win.get(sid, []))
            start = w[0][1] if w else 0
            end = w[-1][1] if w else 0
            ts, te = start >= a.min_repeats, end >= a.min_repeats
            n_both += ts and te
            n_one += ts != te
            fh.write(f"{sid}\t{length}\t{start}\t{end}\t{'yes' if ts else 'no'}\t{'yes' if te else 'no'}\n")
    with open(a.out_summary, "w") as fh:
        fh.write(f"motif\t{a.motif}\n")
        fh.write(f"min_repeats\t{a.min_repeats}\n")
        fh.write(f"sequences\t{len(lengths)}\n")
        fh.write(f"telomere_both_ends\t{n_both}\n")
        fh.write(f"telomere_one_end\t{n_one}\n")


if __name__ == "__main__":
    main()
