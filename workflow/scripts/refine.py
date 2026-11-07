#!/usr/bin/env python3
"""05_EVM/refine: compare EVM models with EviAnn models and build the refined gene set.

Loci = connected groups of EVM and EviAnn coding genes whose CDS overlap on the same strand.
Models are compared by their CDS intron chains (EVM models carry no UTRs, so whole-exon chains
would rarely be identical even for correct models).

  R1  1 EVM : 1 EviAnn, an EviAnn transcript has the same CDS intron chain
        -> keep the EVM CDS, take exons/UTRs from that EviAnn transcript, add other complete
           EviAnn isoforms (refine.import_utr_isoforms)
  R2  EviAnn gene without EVM overlap
        -> add it if it has a complete ORF >= min_cds_aa, RNA support, and CDS-TE overlap
           <= max_te_overlap
  R3  1 EVM : 1 EviAnn with different CDS intron chains
        -> EviAnn chain fully supported by >= min_lr_support long reads and EVM chain not:
           replace EVM with the EviAnn gene; both supported: keep EVM and add the supported EviAnn
           transcripts as isoforms; otherwise keep EVM
  R4  split (1 EVM : n EviAnn) / fusion (n EVM : 1 EviAnn) / n : m
        split : every EviAnn gene supported and no read links them -> replace EVM by the EviAnn
                genes; reads link them -> keep EVM; else keep EVM, needs_review
        fusion: the EviAnn gene is supported -> replace the EVM genes; else keep, needs_review
        n : m : keep EVM, needs_review
  R5  EVM gene without EviAnn overlap
        -> keep if protein (miniprot / GeMoMa), RNA (StringTie) or BRAKER4 gene_support evidence
           overlaps it; otherwise keep but flag refine_support=none (or drop with
           --remove-unsupported-evm)
  R6  EviAnn lncRNAs are split off upstream (split_coding.py) and not merged here

"Supported" with long reads = long reads whose introns inside the model span are exactly the
model's introns (single-exon models: an unspliced read covering >= --single-exon-overlap of the
model). Each read counts with its weight: the full_length_coverage of an Iso-Seq HQ transcript
(number of full-length reads it was built from), 1 for raw reads.

Without long reads (short reads only):
  B  every intron of the model has >= --min-sr-junction junction reads AND at least
     --sr-min-junction-ratio x the reads of the strongest overlapping (competing) intron
  A  a StringTie transcript (short-read assembly, cov >= --sr-stringtie-min-cov, same strand or
     unstranded) whose introns inside the model span are exactly the model's introns
  C  single-exon model: a single-exon StringTie transcript (cov >= --sr-stringtie-min-cov) covers
     >= --single-exon-overlap of it, and no junction with >= --min-sr-junction reads lies inside it
  1 intron: B; >= 2 introns: A and B (A shows the introns belong to one transcript); single exon: C.

Outputs: refined GFF3 (gene attributes refine_rule=, eviann_support=yes|no, refine_support=) and
refine_decisions.tsv (one row per locus).
"""
import argparse
import bisect
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import parse_attrs, read_fasta, read_gff  # noqa: E402
from gff_models import load_models  # noqa: E402

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ---------------------------------------------------------------- helpers kept for other scripts
def cds_len(tx):
    return sum(int(f[4]) - int(f[3]) + 1 for f in tx["feats"] if f[2] == "CDS")


def write_gene(fh, g, prefix="", source=None):
    """Write a load_models() gene unchanged (optionally re-prefixing IDs / changing source)."""
    def fix(row):
        r = row[:]
        if source:
            r[1] = source
        if prefix:
            parts = []
            for kv in r[8].split(";"):
                if kv.startswith(("ID=", "Parent=")):
                    k, v = kv.split("=", 1)
                    kv = k + "=" + ",".join(prefix + x for x in v.split(","))
                parts.append(kv)
            r[8] = ";".join(parts)
        return r
    fh.write("\t".join(fix(g["row"])) + "\n")
    for t in g["tx"].values():
        fh.write("\t".join(fix(t["row"])) + "\n")
        for f in t["feats"]:
            fh.write("\t".join(fix(f)) + "\n")


# ---------------------------------------------------------------- model representation
class Tx:
    __slots__ = ("id", "seqid", "strand", "exons", "cds")

    def __init__(self, tid, seqid, strand, exons, cds):
        self.id, self.seqid, self.strand = tid, seqid, strand
        self.cds = sorted(cds)
        self.exons = sorted(exons) if exons else list(self.cds)

    @property
    def span(self):
        return self.exons[0][0], self.exons[-1][1]

    @property
    def introns(self):
        return tuple((a[1] + 1, b[0] - 1) for a, b in zip(self.exons, self.exons[1:]))

    @property
    def cds_introns(self):
        return tuple((a[1] + 1, b[0] - 1) for a, b in zip(self.cds, self.cds[1:]))

    @property
    def cds_aa(self):
        return sum(e - s + 1 for s, e in self.cds) // 3


class Gene:
    def __init__(self, gid, seqid, strand, txs, origin):
        self.id, self.seqid, self.strand, self.txs, self.origin = gid, seqid, strand, txs, origin

    @property
    def cds_union(self):
        return merge([seg for t in self.txs for seg in t.cds])

    @property
    def span(self):
        return min(t.span[0] for t in self.txs), max(t.span[1] for t in self.txs)


def merge(ivs):
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def overlap_bp(a, b):
    """Overlap of two merged, sorted interval lists."""
    i = j = n = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo <= hi:
            n += hi - lo + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return n


def genes_from(path, origin):
    out = []
    models, _ = load_models(path)
    for gid, g in models.items():
        seqid, strand = g["row"][0], g["row"][6]
        txs = []
        for tid, t in g["tx"].items():
            cds = [(int(f[3]), int(f[4])) for f in t["feats"] if f[2] == "CDS"]
            if not cds:
                continue
            exons = [(int(f[3]), int(f[4])) for f in t["feats"] if f[2] == "exon"]
            txs.append(Tx(tid, seqid, strand, exons, cds))
        if txs:
            out.append(Gene(gid, seqid, strand, txs, origin))
    return out


# ---------------------------------------------------------------- sequence checks
def revcomp(s):
    return s.translate(COMP)[::-1]


def complete_orf(genome, tx):
    seq = genome.get(tx.seqid)
    if seq is None:
        return False
    cds = "".join(seq[s - 1:e] for s, e in tx.cds).upper()
    if tx.strand == "-":
        cds = revcomp(cds)
    if len(cds) % 3 or len(cds) < 6 or cds[:3] != "ATG" or cds[-3:] not in STOPS:
        return False
    return not any(cds[i:i + 3] in STOPS for i in range(0, len(cds) - 3, 3))


# ---------------------------------------------------------------- RNA support
class LongReadSupport:
    kind = "long"
    BIN = 10000

    def __init__(self, path, single_exon_overlap):
        self.reads, self.by_intron, self.bins = [], defaultdict(list), defaultdict(list)
        self.min_ov = single_exon_overlap
        self.weights, self.weighted = [], False
        with open(path) as fh:
            for line in fh:
                if line.startswith("seqid\t"):          # header (repeated when per-sample tables are concatenated)
                    continue
                p = line.rstrip("\n").split("\t")
                seqid, s, e, intr = p[:4]
                w = int(p[4]) if len(p) > 4 and p[4] else 1
                self.weighted |= w != 1
                introns = tuple(tuple(map(int, x.split("-"))) for x in intr.split(",") if x)
                k = len(self.reads)
                self.reads.append((seqid, int(s), int(e), introns))
                self.weights.append(w)
                for it in introns:
                    self.by_intron[(seqid,) + it].append(k)
                for b in range(int(s) // self.BIN, int(e) // self.BIN + 1):
                    self.bins[(seqid, b)].append(k)

    def _near(self, seqid, s, e):
        seen = set()
        for b in range(s // self.BIN, e // self.BIN + 1):
            for k in self.bins.get((seqid, b), ()):
                if k not in seen:
                    seen.add(k)
                    yield k, self.reads[k]

    def support(self, seqid, exons, strand=None):
        exons = sorted(exons)
        s, e = exons[0][0], exons[-1][1]
        chain = tuple((a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:]))
        n = 0
        if chain:
            for k in self.by_intron.get((seqid,) + chain[0], ()):
                rs, rintr = self.reads[k][0], self.reads[k][3]
                inside = tuple(it for it in rintr if it[1] >= s and it[0] <= e)
                if rs == seqid and inside == chain:
                    n += self.weights[k]
            return n
        length = e - s + 1
        for k, (rs, r0, r1, rintr) in self._near(seqid, s, e):
            if any(it[1] >= s and it[0] <= e for it in rintr):
                continue
            if min(e, r1) - max(s, r0) + 1 >= self.min_ov * length:
                n += self.weights[k]
        return n

    def links(self, seqid, a, b):
        """Reads with an intron from inside gene span a to inside gene span b (a before b)."""
        n = 0
        for k, (rs, r0, r1, rintr) in self._near(seqid, a[0], b[1]):
            if any(a[0] <= it[0] - 1 <= a[1] and b[0] <= it[1] + 1 <= b[1] for it in rintr):
                n += self.weights[k]
        return n


LongReadSupport.basis = lambda self, seqid, exons, strand=None: (
    "long" if self.support(seqid, exons, strand) else "-")


class JunctionSupport:
    """Junction read counts (seqid, start, end -> reads); used for links and for criterion B."""
    kind = "short"

    def __init__(self, path, min_reads):
        self.j, self.min = {}, min_reads
        with open(path) as fh:
            next(fh)
            for line in fh:
                seqid, s, e, n = line.rstrip("\n").split("\t")
                self.j[(seqid, int(s), int(e))] = int(n)

    def support(self, seqid, exons, strand=None):
        exons = sorted(exons)
        chain = [(a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:])]
        if not chain:
            return 0
        n = min(self.j.get((seqid,) + it, 0) for it in chain)
        return n if n >= self.min else 0

    def links(self, seqid, a, b):
        return sum(n for (sq, s, e), n in self.j.items()
                   if sq == seqid and a[0] <= s - 1 <= a[1] and b[0] <= e + 1 <= b[1])


class ShortReadSupport(JunctionSupport):
    """Short-read support without long reads: criteria A (StringTie chain), B (junction reads,
    absolute and relative to competing introns) and C (single-exon models); see the module doc."""
    kind = "short (A: StringTie chain, B: junctions, C: single exon)"
    BIN = 10000

    def __init__(self, junctions, stringtie, min_reads, min_ratio, min_cov, single_exon_overlap):
        super().__init__(junctions, min_reads)
        self.ratio, self.min_cov, self.min_ov = min_ratio, min_cov, single_exon_overlap
        self.jbins = defaultdict(list)                       # (seqid, bin) -> [(s, e, n)]
        for (sq, js, je), n in self.j.items():
            for b in range(js // self.BIN, je // self.BIN + 1):
                self.jbins[(sq, b)].append((js, je, n))
        self.tx, self.tbins = [], defaultdict(list)          # StringTie transcripts with cov >= min_cov
        exons, info = defaultdict(list), {}
        if stringtie and os.path.exists(stringtie):
            for c in read_gff(stringtie):
                at = parse_attrs(c[8])
                tid = at.get("transcript_id") or at.get("Parent") or at.get("ID")
                if c[2] == "transcript":
                    try:
                        cov = float(at.get("cov", 0))
                    except ValueError:
                        cov = 0.0
                    info[tid] = (c[0], c[6], cov)
                elif c[2] == "exon" and tid:
                    exons[tid].append((int(c[3]), int(c[4])))
        for tid, ex in exons.items():
            seqid, strand, cov = info.get(tid, (None, ".", 0.0))
            if seqid is None or cov < self.min_cov:
                continue
            ex = sorted(ex)
            introns = tuple((a[1] + 1, b[0] - 1) for a, b in zip(ex, ex[1:]))
            k = len(self.tx)
            self.tx.append((seqid, strand, ex[0][0], ex[-1][1], introns))
            for b in range(ex[0][0] // self.BIN, ex[-1][1] // self.BIN + 1):
                self.tbins[(seqid, b)].append(k)

    def _junctions_near(self, seqid, s, e):
        seen = set()
        for b in range(s // self.BIN, e // self.BIN + 1):
            for x in self.jbins.get((seqid, b), ()):
                if x not in seen:
                    seen.add(x)
                    yield x

    def _tx_near(self, seqid, s, e):
        seen = set()
        for b in range(s // self.BIN, e // self.BIN + 1):
            for k in self.tbins.get((seqid, b), ()):
                if k not in seen:
                    seen.add(k)
                    yield self.tx[k]

    def intron_ok(self, seqid, it):
        """B for one intron: absolute reads and share of the strongest overlapping intron."""
        n = self.j.get((seqid,) + it, 0)
        if n < self.min:
            return False
        top = max((m for js, je, m in self._junctions_near(seqid, it[0], it[1])
                   if js <= it[1] and je >= it[0]), default=n)
        return n >= self.ratio * top

    def chain_in_stringtie(self, seqid, strand, s, e, chain):
        """A: a StringTie transcript whose introns inside [s, e] are exactly the chain."""
        for sq, st, ts, te, introns in self._tx_near(seqid, s, e):
            if sq != seqid or (strand and st not in (strand, ".")):
                continue
            if tuple(it for it in introns if it[1] >= s and it[0] <= e) == chain:
                return True
        return False

    def single_exon_ok(self, seqid, strand, s, e):
        """C: a single-exon StringTie transcript covers the model and no junction lies inside."""
        if any(js > s and je < e and n >= self.min for js, je, n in self._junctions_near(seqid, s, e)):
            return False
        length = e - s + 1
        return any(sq == seqid and not introns and (not strand or st in (strand, "."))
                   and min(e, te) - max(s, ts) + 1 >= self.min_ov * length
                   for sq, st, ts, te, introns in self._tx_near(seqid, s, e))

    def basis(self, seqid, exons, strand=None):
        exons = sorted(exons)
        s, e = exons[0][0], exons[-1][1]
        chain = tuple((a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:]))
        if not chain:
            return "C" if self.single_exon_ok(seqid, strand, s, e) else "-"
        if not all(self.intron_ok(seqid, it) for it in chain):
            return "-"
        if len(chain) == 1:
            return "B"
        return "A+B" if self.chain_in_stringtie(seqid, strand, s, e, chain) else "-"

    def support(self, seqid, exons, strand=None):
        """Weakest intron's junction reads when supported (1 for single-exon models), else 0."""
        b = self.basis(seqid, exons, strand)
        if b == "-":
            return 0
        if b == "C":
            return 1
        exons = sorted(exons)
        return min(self.j.get((seqid, a[1] + 1, c[0] - 1), 0) for a, c in zip(exons, exons[1:]))


class NoSupport:
    kind = "none"

    def support(self, seqid, exons, strand=None):
        return 0

    def basis(self, seqid, exons, strand=None):
        return "-"

    def links(self, seqid, a, b):
        return 0


# ---------------------------------------------------------------- evidence for R2 / R5
class Index:
    """Merged intervals per key for overlap queries."""

    def __init__(self, items):
        g = defaultdict(list)
        for key, s, e in items:
            g[key].append((s, e))
        self.iv = {k: merge(v) for k, v in g.items()}
        self.st = {k: [x[0] for x in v] for k, v in self.iv.items()}

    def covered(self, key, s, e):
        iv, st = self.iv.get(key), self.st.get(key)
        if not iv:
            return 0
        i = max(0, bisect.bisect_right(st, e) - 1)
        n = 0
        while i >= 0 and iv[i][1] >= s:
            n += max(0, min(e, iv[i][1]) - max(s, iv[i][0]) + 1)
            i -= 1
        return n


def gff_index(path, types, stranded=True, skip=()):
    items = []
    if path and os.path.exists(path):
        for c in read_gff(path):
            if c[2] in skip or (types and c[2] not in types):
                continue
            key = (c[0], c[6]) if stranded else c[0]
            items.append((key, int(c[3]), int(c[4])))
    return Index(items)


def braker4_index(gff, support_tsv):
    from filter_final import braker4_supported  # lazy: filter_final imports this module
    ix, note = braker4_supported(support_tsv, gff)
    return ix, note


# ---------------------------------------------------------------- output
def emit(fh, gid, seqid, strand, source, txs, attrs):
    """Write gene / mRNA / exon / CDS (phases recomputed) for a list of (tx_id, exons, cds)."""
    s = min(min(x[1])[0] for x in txs)
    e = max(max(x[1])[1] for x in txs)
    extra = "".join(f";{k}={v}" for k, v in attrs.items())
    fh.write(f"{seqid}\t{source}\tgene\t{s}\t{e}\t.\t{strand}\t.\tID={gid}{extra}\n")
    for tid, exons, cds in txs:
        exons, cds = sorted(exons), sorted(cds)
        fh.write(f"{seqid}\t{source}\tmRNA\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\tID={tid};Parent={gid}\n")
        for i, (a, b) in enumerate(exons, 1):
            fh.write(f"{seqid}\t{source}\texon\t{a}\t{b}\t.\t{strand}\t.\tID={tid}.exon{i};Parent={tid}\n")
        order = cds if strand != "-" else list(reversed(cds))
        done, phases = 0, {}
        for a, b in order:
            phases[(a, b)] = (3 - done % 3) % 3
            done += b - a + 1
        for a, b in cds:
            fh.write(f"{seqid}\t{source}\tCDS\t{a}\t{b}\t.\t{strand}\t{phases[(a, b)]}\tID={tid}.cds;Parent={tid}\n")


# ---------------------------------------------------------------- loci
def build_loci(evm, evi):
    allg = evm + evi
    parent = list(range(len(allg)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    order = sorted(range(len(allg)), key=lambda i: (allg[i].seqid, allg[i].strand, allg[i].span[0]))
    active = []
    unions = [g.cds_union for g in allg]
    for i in order:
        g = allg[i]
        active = [j for j in active if allg[j].seqid == g.seqid and allg[j].strand == g.strand
                  and allg[j].span[1] >= g.span[0]]
        for j in active:
            if allg[j].origin != g.origin and overlap_bp(unions[i], unions[j]) > 0:
                parent[find(i)] = find(j)
        active.append(i)
    groups = defaultdict(list)
    for i in range(len(allg)):
        groups[find(i)].append(allg[i])
    return sorted(groups.values(), key=lambda l: (l[0].seqid, min(g.span[0] for g in l)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--evm", required=True)
    ap.add_argument("--eviann", required=True, help="EviAnn coding genes")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--lr-chains", default="", help="extract_introns.py chains output (long reads)")
    ap.add_argument("--sr-junctions", default="", help="extract_introns.py junctions output (fallback)")
    ap.add_argument("--stringtie", default="")
    ap.add_argument("--proteins", nargs="*", default=[], help="miniprot / GeMoMa GFF3 (CDS used)")
    ap.add_argument("--braker4-gff", default="")
    ap.add_argument("--braker4-support", default="")
    ap.add_argument("--te-gff", default="")
    ap.add_argument("--min-cds-aa", type=int, default=50)
    ap.add_argument("--max-te-overlap", type=float, default=0.5)
    ap.add_argument("--min-lr-support", type=int, default=2)
    ap.add_argument("--min-sr-junction", type=int, default=3)
    ap.add_argument("--sr-min-junction-ratio", type=float, default=0.1)
    ap.add_argument("--sr-stringtie-min-cov", type=float, default=3.0)
    ap.add_argument("--single-exon-overlap", type=float, default=0.5)
    ap.add_argument("--import-utr-isoforms", choices=["true", "false"], default="true")
    ap.add_argument("--remove-unsupported-evm", choices=["true", "false"], default="false")
    ap.add_argument("--out", required=True)
    ap.add_argument("--decisions", required=True)
    a = ap.parse_args()

    genome = {h.split()[0]: s for h, s in read_fasta(a.genome)}
    if a.lr_chains and os.path.exists(a.lr_chains):
        sup = LongReadSupport(a.lr_chains, a.single_exon_overlap)
        min_sup = a.min_lr_support
    elif a.sr_junctions and os.path.exists(a.sr_junctions):
        sup = ShortReadSupport(a.sr_junctions, a.stringtie, a.min_sr_junction, a.sr_min_junction_ratio,
                               a.sr_stringtie_min_cov, a.single_exon_overlap)
        min_sup = 1
    else:
        sup, min_sup = NoSupport(), 1
    st_ix = gff_index(a.stringtie, {"exon"})
    prot_ix = Index([((c[0], c[6]), int(c[3]), int(c[4])) for p in a.proteins if os.path.exists(p)
                     for c in read_gff(p) if c[2] == "CDS"])
    te_ix = gff_index(a.te_gff, None, stranded=False,
                      skip={"repeat_region", "long_terminal_repeat", "target_site_duplication"})
    b4_ix, b4_note = braker4_index(a.braker4_gff, a.braker4_support)
    imp = a.import_utr_isoforms == "true"

    evm, evi = genes_from(a.evm, "EVM"), genes_from(a.eviann, "EviAnn")
    loci = build_loci(evm, evi)

    def supported(tx):
        return sup.support(tx.seqid, tx.exons, tx.strand)

    def basis(genes):
        b = sorted({sup.basis(t.seqid, t.exons, t.strand) for g in genes for t in g.txs} - {"-"})
        return "/".join(b) or "-"

    def good_eviann_txs(g):
        return [t for t in g.txs if t.cds_aa >= a.min_cds_aa and complete_orf(genome, t)]

    def rna_overlap(g):
        return any(st_ix.covered((g.seqid, g.strand), s, e) or st_ix.covered((g.seqid, "."), s, e)
                   for t in g.txs for s, e in t.exons)

    def te_frac(g):
        u = g.cds_union
        tot = sum(e - s + 1 for s, e in u)
        return sum(te_ix.covered(g.seqid, s, e) for s, e in u) / tot if tot else 0.0

    out = open(a.out, "w")
    dec = open(a.decisions, "w")
    out.write("##gff-version 3\n")
    weighting = (" (weighted by full_length_coverage)" if getattr(sup, "weighted", False) else "")
    dec.write(f"# support={sup.kind}{weighting} min_support={min_sup} min_cds_aa={a.min_cds_aa} "
              f"max_te_overlap={a.max_te_overlap}; {b4_note}\n")
    dec.write("locus\tseqid\tstart\tend\tstrand\tevm_genes\teviann_genes\trule\taction\t"
              "evm_support\teviann_support\tevidence\tnote\tsupport_basis\n")

    def keep_evm(g, rule, extra=None, eviann="no", support=None):
        attrs = {"refine_rule": rule, "eviann_support": eviann}
        if support:
            attrs["refine_support"] = support
        txs = [(t.id, t.exons, t.cds) for t in g.txs] + (extra or [])
        emit(out, g.id, g.seqid, g.strand, "EVM", txs, attrs)

    def add_eviann(g, rule, txs):
        emit(out, "EviAnn." + g.id, g.seqid, g.strand, "EviAnn_refine",
             [("EviAnn." + t.id, t.exons, t.cds) for t in txs], {"refine_rule": rule, "eviann_support": "yes"})

    for n, locus in enumerate(loci, 1):
        E = [g for g in locus if g.origin == "EVM"]
        V = [g for g in locus if g.origin == "EviAnn"]
        g0 = locus[0]
        s = min(g.span[0] for g in locus)
        e = max(g.span[1] for g in locus)
        evm_sup = max((supported(t) for g in E for t in g.txs), default=0)
        evi_sup = max((supported(t) for g in V for t in g.txs), default=0)
        rule, action, evidence, note = "", "", "-", ""

        if not V:                                             # R5
            g = E[0]
            ev = []
            if any(prot_ix.covered((g.seqid, g.strand), cs, ce) for cs, ce in g.cds_union):
                ev.append("protein")
            if rna_overlap(g) or evm_sup >= min_sup:
                ev.append("RNA")
            if any(b4_ix.covered((g.seqid, g.strand), cs, ce) for cs, ce in g.cds_union):
                ev.append("BRAKER4_support")
            rule, evidence = "R5", ",".join(ev) or "none"
            if ev:
                action = "keep_evm"
                keep_evm(g, rule, support=evidence)
            elif a.remove_unsupported_evm == "true":
                action, note = "remove", "no protein/RNA/BRAKER4 evidence"
            else:
                action, note = "keep_evm_flagged", "no protein/RNA/BRAKER4 evidence"
                keep_evm(g, rule, support="none")

        elif not E:                                           # R2
            g = V[0]
            good = good_eviann_txs(g)
            rna = rna_overlap(g) or evi_sup >= min_sup
            tf = te_frac(g)
            rule, evidence = "R2", f"RNA={'yes' if rna else 'no'};TE={tf:.2f}"
            if not good:
                action, note = "skip", f"no complete ORF >= {a.min_cds_aa} aa"
            elif not rna:
                action, note = "skip", "no RNA support"
            elif tf > a.max_te_overlap:
                action, note = "skip", f"CDS-TE overlap {tf:.2f} > {a.max_te_overlap}"
            else:
                action = "add_eviann"
                add_eviann(g, rule, good)

        elif len(E) == 1 and len(V) == 1:                     # R1 / R3
            ge, gv = E[0], V[0]
            et = ge.txs[0]
            match = [t for t in gv.txs if t.cds_introns == et.cds_introns]
            if match:
                rule = "R1"
                best = match[0]
                compatible = all(any(xs <= cs and ce <= xe for xs, xe in best.exons) for cs, ce in et.cds)
                extra = []
                if imp:
                    used = {(tuple(t.exons), tuple(t.cds)) for t in ge.txs}
                    for k, t in enumerate(good_eviann_txs(gv), 1):
                        if t is best or (tuple(t.exons), tuple(t.cds)) in used:
                            continue
                        extra.append((f"{ge.id}.EviAnn{k}", t.exons, t.cds))
                if compatible and imp:
                    action = "merge_utr" + (f"+{len(extra)}_isoforms" if extra else "")
                    txs = [(et.id, best.exons, et.cds)] + [(t.id, t.exons, t.cds) for t in ge.txs[1:]] + extra
                    emit(out, ge.id, ge.seqid, ge.strand, "EVM", txs, {"refine_rule": rule, "eviann_support": "yes"})
                else:
                    action = "keep_evm"
                    note = "EVM CDS not inside EviAnn exons" if not compatible else "isoform import disabled"
                    keep_evm(ge, rule, eviann="yes")
            else:
                rule = "R3"
                good = good_eviann_txs(gv)
                sup_txs = [t for t in good if supported(t) >= min_sup]
                if sup_txs and evm_sup < min_sup:
                    action = "replace_with_eviann"
                    add_eviann(gv, rule, good)
                elif sup_txs:
                    action = "keep_evm+eviann_isoforms"
                    keep_evm(ge, rule, eviann="yes",
                             extra=[(f"{ge.id}.EviAnn{k}", t.exons, t.cds) for k, t in enumerate(sup_txs, 1)])
                else:
                    action = "keep_evm"
                    note = "EviAnn chain not supported" if good else "no complete EviAnn ORF"
                    keep_evm(ge, rule)

        elif len(E) == 1:                                     # R4 split
            rule = "R4_split"
            ge = E[0]
            goods = [good_eviann_txs(g) for g in V]
            each = all(any(supported(t) >= min_sup for t in gd) for gd in goods)
            vs = sorted(V, key=lambda g: g.span[0])
            linked = sum(sup.links(ge.seqid, x.span, y.span) for x, y in zip(vs, vs[1:]))
            evidence = f"linking_reads={linked}"
            if each and linked == 0:
                action = "replace_with_eviann"
                for g, gd in zip(V, goods):
                    add_eviann(g, rule, gd)
            elif linked > 0:
                action, note = "keep_evm", "reads link the EviAnn genes (fusion supported)"
                keep_evm(ge, rule)
            else:
                action, note = "keep_evm", "needs_review: split not supported"
                keep_evm(ge, rule, support="needs_review")

        elif len(V) == 1:                                     # R4 fusion
            rule = "R4_fusion"
            gv = V[0]
            good = good_eviann_txs(gv)
            spanning = [t for t in good if supported(t) >= min_sup
                        and all(overlap_bp(merge(t.cds), g.cds_union) > 0 for g in E)]
            if spanning:
                action = "replace_with_eviann"
                add_eviann(gv, rule, good)
            else:
                action, note = "keep_evm", "needs_review: fusion not supported"
                for g in E:
                    keep_evm(g, rule, support="needs_review")

        else:                                                 # R4 many : many
            rule, action, note = "R4_complex", "keep_evm", "needs_review"
            for g in E:
                keep_evm(g, rule, support="needs_review")

        dec.write("\t".join(map(str, [
            f"L{n:06d}", g0.seqid, s, e, g0.strand, ",".join(g.id for g in E) or "-",
            ",".join(g.id for g in V) or "-", rule, action, evm_sup, evi_sup, evidence, note or "-",
            f"evm:{basis(E)};eviann:{basis(V)}"])) + "\n")
    out.close()
    dec.close()


if __name__ == "__main__":
    main()
