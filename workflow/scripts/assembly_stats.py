#!/usr/bin/env python3
"""Assembly contiguity statistics at scaffold and contig level (standard library only).

Contigs are obtained by splitting scaffolds at runs of >= --gap-min-n 'N'/'n'. Shorter N runs
stay inside contigs (counted as ambiguous bases, not gaps).
Outputs:
  --out-summary  level-wise metrics (count, total, min/max, N50/L50, N90/L90, GC, N bases, gaps)
  --out-seqs     per scaffold: length, GC, N bases, gaps, contigs (used for plots/telomeres)
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_fasta  # noqa: E402


def nx(lengths, frac):
    """(Nx, Lx) for sorted-descending lengths."""
    total = sum(lengths)
    acc = 0
    for i, l in enumerate(lengths, 1):
        acc += l
        if acc >= total * frac:
            return l, i
    return 0, 0


def level_stats(lengths, gc, n_bases, at_gc_total):
    lengths = sorted(lengths, reverse=True)
    n50, l50 = nx(lengths, 0.5)
    n90, l90 = nx(lengths, 0.9)
    return {
        "count": len(lengths), "total_length": sum(lengths),
        "max_length": lengths[0] if lengths else 0, "min_length": lengths[-1] if lengths else 0,
        "N50": n50, "L50": l50, "N90": n90, "L90": l90,
        "GC_percent": f"{100 * gc / at_gc_total:.2f}" if at_gc_total else "0",
        "N_bases": n_bases,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gap-min-n", type=int, default=10)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-seqs", required=True)
    a = ap.parse_args()
    gap_re = re.compile("[Nn]{%d,}" % a.gap_min_n)

    scaf_len, contig_len, per_seq = [], [], []
    gc = acgt = n_total = gaps = gap_bp = 0
    for header, seq in read_fasta(a.fasta):
        sid = header.split()[0]
        u = seq.upper()
        s_gc = u.count("G") + u.count("C")
        s_at = u.count("A") + u.count("T")
        s_n = u.count("N")
        s_gaps = [m.end() - m.start() for m in gap_re.finditer(seq)]
        pieces = [p for p in gap_re.split(seq) if p]
        scaf_len.append(len(seq))
        contig_len += [len(p) for p in pieces]
        gc += s_gc
        acgt += s_gc + s_at
        n_total += s_n
        gaps += len(s_gaps)
        gap_bp += sum(s_gaps)
        per_seq.append((sid, len(seq), f"{100 * s_gc / (s_gc + s_at):.2f}" if s_gc + s_at else "0",
                        s_n, len(s_gaps), sum(s_gaps), len(pieces)))
    if not scaf_len:
        sys.exit("ERROR: empty FASTA")

    scaf = level_stats(scaf_len, gc, n_total, acgt)
    # contig level: gap Ns removed, short N runs remain inside contigs
    contig = level_stats(contig_len, gc, n_total - gap_bp, acgt)
    with open(a.out_summary, "w") as fh:
        fh.write("metric\tscaffold\tcontig\n")
        for k in scaf:
            fh.write(f"{k}\t{scaf[k]}\t{contig[k]}\n")
        fh.write(f"gaps\t{gaps}\t-\n")
        fh.write(f"gap_length\t{gap_bp}\t-\n")
        fh.write(f"gap_min_n\t{a.gap_min_n}\t-\n")
    with open(a.out_seqs, "w") as fh:
        fh.write("seq_id\tlength\tGC_percent\tN_bases\tgaps\tgap_length\tcontigs\n")
        for r in sorted(per_seq, key=lambda x: -x[1]):
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
