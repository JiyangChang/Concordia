#!/usr/bin/env python3
"""Read SAM records from stdin (e.g. `samtools view -F 0x904 x.bam | extract_introns.py ...`).

  chains     one line per spliced or unspliced read: seqid, start, end, introns ("s-e,s-e,..."),
             weight. Used for long-read support of whole intron chains (refine R3/R4).
             With --reads, the weight of a read is the full_length_coverage=N value in its
             FASTA/FASTQ header (Iso-Seq cluster / polish HQ transcripts: number of full-length
             reads behind the consensus); reads without that field (raw FLNC / ONT reads) weigh 1.
             With --counts, weights come from an Iso-Seq count table instead (for HQ transcripts whose
             headers lost the tag): isoseq cluster's cluster_report.csv (cluster_id,read_id,read_type:
             rows per cluster are counted), or a table whose first column is the transcript ID and
             whose other numeric columns are full-length read counts (summed), e.g. isoseq collapse's
             *.flnc_count.txt or a simple "id<TAB>count" file. CSV or TSV, header optional.
             --summary writes reads / weighted reads / weight source for the report.
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


def read_counts(path):
    """transcript ID -> full-length read count from an Iso-Seq count table (see module doc)."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]
    if not lines:
        return {}
    sep = "," if "," in lines[0] else None
    rows = [l.split(sep) for l in lines]
    head = [x.strip().lower() for x in rows[0]]
    if "cluster_id" in head and "read_id" in head:             # isoseq cluster: one row per FL read
        k = head.index("cluster_id")
        return dict(Counter(r[k].strip() for r in rows[1:] if len(r) > k))

    def num(x):
        try:
            return float(x)
        except ValueError:
            return None
    body = rows[1:] if all(num(x) is None for x in rows[0][1:]) else rows   # skip a header line
    out = {}
    for r in body:
        vals = [v for v in (num(x) for x in r[1:]) if v is not None]
        if vals:
            out[r[0].strip()] = int(round(sum(vals)))
    return out


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
    ap.add_argument("--counts", default="", help="Iso-Seq count table (chains: read weights; overrides headers)")
    ap.add_argument("--summary", default="", help="chains: write reads / weighted reads / weight source")
    ap.add_argument("--sample", default="", help="sample name for --summary")
    ap.add_argument("--platform", default="", help="platform for --summary")
    a = ap.parse_args()
    counts = Counter()
    weights, source = {}, "none"
    if a.mode == "chains" and a.reads:
        weights = read_weights(a.reads)
        source = "full_length_coverage in read headers"
        print(f"{len(weights)} reads carry full_length_coverage", file=sys.stderr)
    if a.mode == "chains" and a.counts:
        weights = read_counts(a.counts)
        source = f"count table {a.counts.rsplit('/', 1)[-1]}"
        print(f"{len(weights)} transcripts with full-length read counts in {a.counts}", file=sys.stderr)
    n_reads = n_weighted = total_w = 0
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
                w = weights.get(name)
                n_reads += 1
                n_weighted += w is not None
                total_w += w if w is not None else 1
                out.write(f"{seqid}\t{pos}\t{end}\t{','.join(f'{s}-{e}' for s, e in introns)}\t{w or 1}\n")
            else:
                for s, e in introns:
                    counts[(seqid, s, e)] += 1
        if a.mode == "junctions":
            out.write("seqid\tstart\tend\treads\n")
            for (seqid, s, e), n in sorted(counts.items()):
                out.write(f"{seqid}\t{s}\t{e}\t{n}\n")
    if a.mode == "chains" and a.summary:
        with open(a.summary, "w") as fh:
            fh.write("sample\tplatform\taligned_reads\tweighted_reads\tweighted_%\tfull_length_reads\tweight_source\n")
            pct = f"{100 * n_weighted / n_reads:.1f}" if n_reads else "0.0"
            fh.write(f"{a.sample}\t{a.platform}\t{n_reads}\t{n_weighted}\t{pct}\t{total_w}\t{source}\n")


if __name__ == "__main__":
    main()
