#!/usr/bin/env python3
"""Read SAM records from stdin (e.g. `samtools view -F 0x904 x.bam | extract_introns.py ...`).

  chains     one line per spliced or unspliced read: seqid, start, end, introns ("s-e,s-e,..."),
             weight. Used for long-read support of whole intron chains (refine R3/R4).
             With --reads, the weight of a read is the full_length_coverage=N value in its
             FASTA/FASTQ header (Iso-Seq cluster / polish HQ transcripts: number of full-length
             reads behind the consensus); reads without that field (raw FLNC / ONT reads) weigh 1.
  junctions  aggregated intron counts: seqid, start, end, reads (short-read fallback)
Coordinates are 1-based, inclusive; introns come from CIGAR 'N' operations.
"""
import argparse
import gzip
import re
import sys
from collections import Counter

FLC = re.compile(r"full_length_coverage=(\d+)")


def read_weights(path):
    """read name -> full_length_coverage, for headers that carry it (FASTA or FASTQ, gz ok)."""
    w = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        fastq = None
        for i, line in enumerate(fh):
            if fastq is None:
                fastq = line.startswith("@")
            if (fastq and i % 4 == 0) or (not fastq and line.startswith(">")):
                m = FLC.search(line)
                if m:
                    w[line[1:].split()[0]] = int(m.group(1))
    return w


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = set("MDN=X")


def blocks(pos, cigar):
    """(alignment end, list of introns) for a read starting at 1-based pos."""
    introns, p = [], pos
    for n, op in CIGAR.findall(cigar):
        n = int(n)
        if op == "N":
            introns.append((p, p + n - 1))
        if op in REF_CONSUMING:
            p += n
    return p - 1, introns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["chains", "junctions"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--reads", default="", help="FASTA/FASTQ the BAM was made from (chains: read weights)")
    a = ap.parse_args()
    counts = Counter()
    weights = read_weights(a.reads) if a.mode == "chains" and a.reads else {}
    if a.mode == "chains" and a.reads:
        print(f"{len(weights)} reads carry full_length_coverage", file=sys.stderr)
    with open(a.out, "w") as out:
        if a.mode == "chains":
            out.write("seqid\tstart\tend\tintrons\tweight\n")
        for line in sys.stdin:
            if line.startswith("@"):
                continue
            f = line.split("\t", 6)
            if len(f) < 6 or f[5] == "*":
                continue
            name, seqid, pos, cigar = f[0], f[2], int(f[3]), f[5]
            end, introns = blocks(pos, cigar)
            if a.mode == "chains":
                out.write(f"{seqid}\t{pos}\t{end}\t{','.join(f'{s}-{e}' for s, e in introns)}\t{weights.get(name, 1)}\n")
            else:
                for s, e in introns:
                    counts[(seqid, s, e)] += 1
        if a.mode == "junctions":
            out.write("seqid\tstart\tend\treads\n")
            for (seqid, s, e), n in sorted(counts.items()):
                out.write(f"{seqid}\t{s}\t{e}\t{n}\n")


if __name__ == "__main__":
    main()
