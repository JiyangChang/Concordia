"""refine.py: rules R1-R5 on a synthetic genome with real start / stop codons."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "refine.py")
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
L = 30000


def put_orf(seq, cds, strand, start="ATG", stop="TAA"):
    """Write a start codon at the CDS start and a stop codon at its end (filler 'C' has no stops)."""
    pos = [p for s, e in sorted(cds) for p in range(s, e + 1)]
    if strand == "-":
        pos = pos[::-1]
    for i, base in list(enumerate(start)) + [(len(pos) - 3 + i, b) for i, b in enumerate(stop)]:
        seq[pos[i] - 1] = base if strand == "+" else COMP[base]


def gene(gid, cds, exons=None, strand="+", src="x", txs=None):
    """txs: list of (tid, exons, cds) for multi-isoform genes."""
    txs = txs or [(gid + ".t1", exons or cds, cds)]
    lo = min(x[0] for t in txs for x in t[1])
    hi = max(x[1] for t in txs for x in t[1])
    out = [f"c1\t{src}\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID={gid}"]
    for tid, ex, cd in txs:
        out.append(f"c1\t{src}\tmRNA\t{ex[0][0]}\t{ex[-1][1]}\t.\t{strand}\t.\tID={tid};Parent={gid}")
        out += [f"c1\t{src}\texon\t{s}\t{e}\t.\t{strand}\t.\tParent={tid}" for s, e in ex]
        out += [f"c1\t{src}\tCDS\t{s}\t{e}\t.\t{strand}\t0\tParent={tid}" for s, e in cd]
    return "\n".join(out) + "\n"


def read_line(start, end, introns=(), weight=None):
    w = f"\t{weight}" if weight is not None else ""
    return f"c1\t{start}\t{end}\t{','.join(f'{s}-{e}' for s, e in introns)}{w}\n"


def run(tmp_path, evm, evi, orfs, reads=(), junctions=None, stringtie=(), miniprot=(), te=(), extra=()):
    seq = ["C"] * L
    for cds, strand in orfs:
        put_orf(seq, cds, strand)
    p = lambda n: str(tmp_path / n)
    open(p("g.fa"), "w").write(">c1\n" + "".join(seq) + "\n")
    open(p("evm.gff3"), "w").write("".join(evm))
    open(p("evi.gff3"), "w").write("".join(evi))
    open(p("st.gtf"), "w").write("".join(f'c1\tStringTie\texon\t{x[0]}\t{x[1]}\t.\t{x[2] if len(x) > 2 else "+"}\t.\t'
                                         f'transcript_id "S{i}";\n' for i, x in enumerate(stringtie)))
    open(p("mp.gff3"), "w").write("".join(f"c1\tminiprot\tCDS\t{s}\t{e}\t.\t+\t0\tParent=MP{i}\n"
                                          for i, (s, e) in enumerate(miniprot)))
    open(p("te.gff3"), "w").write("".join(f"c1\tEDTA\tLTR_retrotransposon\t{s}\t{e}\t.\t+\t.\tID=te{i}\n"
                                          for i, (s, e) in enumerate(te)))
    args = [sys.executable, SCRIPT, "--evm", p("evm.gff3"), "--eviann", p("evi.gff3"), "--genome", p("g.fa"),
            "--stringtie", p("st.gtf"), "--proteins", p("mp.gff3"), "--te-gff", p("te.gff3"),
            "--out", p("out.gff3"), "--decisions", p("dec.tsv"), *extra]
    if junctions is not None:
        open(p("j.tsv"), "w").write("seqid\tstart\tend\treads\n" +
                                    "".join(f"c1\t{s}\t{e}\t{n}\n" for s, e, n in junctions))
        args += ["--sr-junctions", p("j.tsv")]
    else:
        open(p("lr.tsv"), "w").write("seqid\tstart\tend\tintrons\n" + "".join(reads))
        args += ["--lr-chains", p("lr.tsv")]
    r = subprocess.run(args, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    dec = list(csv.DictReader((l for l in open(p("dec.tsv")) if not l.startswith("#")), delimiter="\t"))
    rows = [l.rstrip("\n").split("\t") for l in open(p("out.gff3")) if not l.startswith("#")]
    return dec, rows


def genes_out(rows):
    return {dict(kv.split("=", 1) for kv in r[8].split(";"))["ID"]: r for r in rows if r[2] == "gene"}


def attrs(row):
    return dict(kv.split("=", 1) for kv in row[8].split(";"))


# ------------------------------------------------------------------------------------------- R1
def test_r1_merges_utr_and_imports_isoform(tmp_path):
    cds = [(1001, 1300), (1501, 1800)]
    evm = [gene("evm1", cds, src="EVM")]
    evi = [gene("v1", None, txs=[("v1.t1", [(901, 1300), (1501, 1900)], cds),
                                 ("v1.t2", [(951, 1300), (1501, 1850)], cds)])]
    dec, rows = run(tmp_path, evm, evi, [(cds, "+")])
    assert dec[0]["rule"] == "R1" and dec[0]["action"] == "merge_utr+1_isoforms"
    g = genes_out(rows)["evm1"]
    assert attrs(g)["eviann_support"] == "yes" and (g[3], g[4]) == ("901", "1900")
    mrnas = [r for r in rows if r[2] == "mRNA"]
    assert len(mrnas) == 2 and (mrnas[0][3], mrnas[0][4]) == ("901", "1900")
    assert [(r[3], r[4]) for r in rows if r[2] == "CDS" and "evm1.t1" in r[8]] == [("1001", "1300"), ("1501", "1800")]


# ------------------------------------------------------------------------------------------- R2
def test_r2_add_and_skip_reasons(tmp_path):
    ok, no_rna, te_hit, short, broken = [(3001, 3300)], [(4001, 4300)], [(5001, 5300)], [(6001, 6090)], [(7001, 7300)]
    evi = [gene("add", ok), gene("norna", no_rna), gene("te", te_hit), gene("short", short), gene("broken", broken)]
    orfs = [(ok, "+"), (no_rna, "+"), (te_hit, "+"), (short, "+")]          # 'broken' has no start codon
    dec, rows = run(tmp_path, [], evi, orfs, stringtie=[(3000, 3400), (5000, 5400), (6000, 6100), (7000, 7400)],
                    te=[(5001, 5250)])
    d = {x["eviann_genes"]: x for x in dec}
    assert d["add"]["action"] == "add_eviann"
    assert d["norna"]["note"] == "no RNA support"
    assert d["te"]["note"].startswith("CDS-TE overlap 0.83")
    assert d["short"]["note"].startswith("no complete ORF") and d["broken"]["note"].startswith("no complete ORF")
    assert list(genes_out(rows)) == ["EviAnn.add"]


def test_r2_minus_strand_orf_and_phases(tmp_path):
    cds = [(3001, 3100), (3201, 3400)]                                       # 100 + 200 bp
    dec, rows = run(tmp_path, [], [gene("m", cds, strand="-")], [(cds, "-")], stringtie=[(3000, 3500, "-")])
    assert dec[0]["action"] == "add_eviann"
    ph = {(r[3], r[4]): r[7] for r in rows if r[2] == "CDS"}
    assert ph[("3201", "3400")] == "0" and ph[("3001", "3100")] == "1"       # 200 bp first -> next phase 1


# ------------------------------------------------------------------------------------------- R3
EVM3 = [(5001, 5300), (5501, 5800)]
EVI3 = [(5001, 5300), (5601, 5900)]


def test_r3_replace_when_only_eviann_supported(tmp_path):
    reads = [read_line(4950, 5950, [(5301, 5600)])] * 2
    dec, rows = run(tmp_path, [gene("e3", EVM3, src="EVM")], [gene("v3", EVI3)], [(EVI3, "+"), (EVM3, "+")], reads)
    assert dec[0]["rule"] == "R3" and dec[0]["action"] == "replace_with_eviann"
    assert dec[0]["eviann_support"] == "2" and dec[0]["evm_support"] == "0"
    assert list(genes_out(rows)) == ["EviAnn.v3"]


def test_r3_keep_without_support(tmp_path):
    dec, rows = run(tmp_path, [gene("e3", EVM3, src="EVM")], [gene("v3", EVI3)], [(EVI3, "+")])
    assert dec[0]["action"] == "keep_evm" and list(genes_out(rows)) == ["e3"]


def test_r3_both_supported_adds_isoform(tmp_path):
    reads = [read_line(4950, 5950, [(5301, 5600)])] * 2 + [read_line(4950, 5850, [(5301, 5500)])] * 2
    dec, rows = run(tmp_path, [gene("e3", EVM3, src="EVM")], [gene("v3", EVI3)], [(EVI3, "+"), (EVM3, "+")], reads)
    assert dec[0]["action"] == "keep_evm+eviann_isoforms"
    assert len([r for r in rows if r[2] == "mRNA"]) == 2


def test_r3_read_with_extra_intron_does_not_support(tmp_path):
    # the read has an additional intron inside the EviAnn model span -> not the same chain
    reads = [read_line(4950, 5950, [(5101, 5150), (5301, 5600)])] * 3
    dec, _ = run(tmp_path, [gene("e3", EVM3, src="EVM")], [gene("v3", EVI3)], [(EVI3, "+")], reads)
    assert dec[0]["eviann_support"] == "0" and dec[0]["action"] == "keep_evm"


def test_r3_short_read_fallback(tmp_path):
    dec, _ = run(tmp_path, [gene("e3", EVM3, src="EVM")], [gene("v3", EVI3)], [(EVI3, "+")],
                 junctions=[(5301, 5600, 10)])
    assert dec[0]["action"] == "replace_with_eviann"
    assert dec[0]["support_basis"] == "evm:-;eviann:B"            # one intron -> criterion B
    dec, _ = run(tmp_path, [gene("e3", EVM3, src="EVM")], [gene("v3", EVI3)], [(EVI3, "+")],
                 junctions=[(5301, 5600, 2)])                                 # < 3 reads
    assert dec[0]["action"] == "keep_evm"


# ------------------------------------------------------------------------------------------- R4
A, B = [(7001, 7300)], [(7701, 8000)]
A_EX, B_EX = [(6951, 7350)], [(7651, 8050)]


def split_case(tmp_path, reads):
    evm = [gene("e4", [(7001, 7300), (7701, 8000)], src="EVM")]
    evi = [gene("vA", A, A_EX), gene("vB", B, B_EX)]
    return run(tmp_path, evm, evi, [(A, "+"), (B, "+")], reads)


def test_r4_split_replaces_when_supported_and_unlinked(tmp_path):
    reads = [read_line(6960, 7340)] * 2 + [read_line(7660, 8040)] * 2
    dec, rows = split_case(tmp_path, reads)
    assert dec[0]["rule"] == "R4_split" and dec[0]["action"] == "replace_with_eviann"
    assert sorted(genes_out(rows)) == ["EviAnn.vA", "EviAnn.vB"]


def test_r4_split_kept_when_reads_link_genes(tmp_path):
    reads = [read_line(6960, 7340)] * 2 + [read_line(7660, 8040)] * 2 + [read_line(7000, 8000, [(7351, 7650)])]
    dec, rows = split_case(tmp_path, reads)
    assert dec[0]["action"] == "keep_evm" and "fusion supported" in dec[0]["note"]


def test_r4_split_without_support_needs_review(tmp_path):
    dec, rows = split_case(tmp_path, [])
    assert "needs_review" in dec[0]["note"] and attrs(genes_out(rows)["e4"])["refine_support"] == "needs_review"


def test_r4_fusion(tmp_path):
    cds = [(9001, 9300), (9601, 9900)]
    evm = [gene("eA", [(9001, 9300)], src="EVM"), gene("eB", [(9601, 9900)], src="EVM")]
    reads = [read_line(8950, 9950, [(9301, 9600)])] * 2
    dec, rows = run(tmp_path, evm, [gene("vF", cds)], [(cds, "+")], reads)
    assert dec[0]["rule"] == "R4_fusion" and dec[0]["action"] == "replace_with_eviann"
    assert list(genes_out(rows)) == ["EviAnn.vF"]
    dec, rows = run(tmp_path, evm, [gene("vF", cds)], [(cds, "+")], [])
    assert dec[0]["action"] == "keep_evm" and sorted(genes_out(rows)) == ["eA", "eB"]


# ------------------------------------------------------------------------------------------- R5
def test_r5_evidence_and_flag(tmp_path):
    evm = [gene("prot", [(11001, 11300)], src="EVM"), gene("bare", [(12001, 12300)], src="EVM")]
    dec, rows = run(tmp_path, evm, [], [], miniprot=[(11050, 11200)])
    d = {x["evm_genes"]: x for x in dec}
    assert d["prot"]["action"] == "keep_evm" and d["prot"]["evidence"] == "protein"
    assert d["bare"]["action"] == "keep_evm_flagged"
    assert attrs(genes_out(rows)["bare"])["refine_support"] == "none"
    dec, rows = run(tmp_path, evm, [], [], miniprot=[(11050, 11200)], extra=["--remove-unsupported-evm", "true"])
    assert list(genes_out(rows)) == ["prot"]


def test_opposite_strand_not_same_locus(tmp_path):
    cds = [(13001, 13300)]
    dec, _ = run(tmp_path, [gene("eP", cds, src="EVM")], [gene("vM", cds, strand="-")], [(cds, "-")],
                 stringtie=[(13000, 13400, "-")])
    assert sorted(x["rule"] for x in dec) == ["R2", "R5"]


# ------------------------------------------------------------------------------------ extract_introns
def test_extract_introns(tmp_path):
    ex = os.path.join(ROOT, "workflow", "scripts", "extract_introns.py")
    sam = ("@HD\tVN:1.6\n"
           "r1\t0\tc1\t100\t60\t10S50M200N30M2D20M100N40M5S\t*\t0\t0\tACGT\t*\n"   # introns 150-349, 402-501
           "r2\t0\tc1\t100\t60\t50M200N50M\t*\t0\t0\tACGT\t*\n"
           "r3\t0\tc1\t5000\t60\t100M\t*\t0\t0\tACGT\t*\n")
    run_ = lambda mode, out: subprocess.run([sys.executable, ex, mode, "--out", str(out)], input=sam,
                                            capture_output=True, text=True, check=True)
    run_("chains", tmp_path / "c.tsv")
    lines = (tmp_path / "c.tsv").read_text().splitlines()
    assert lines[1] == "c1\t100\t541\t150-349,402-501\t1"     # 50M 200N 30M 2D 20M 100N 40M -> ends at 541
    assert lines[3] == "c1\t5000\t5099\t\t1"
    run_("junctions", tmp_path / "j.tsv")
    assert (tmp_path / "j.tsv").read_text().splitlines()[1:] == ["c1\t150\t349\t2", "c1\t402\t501\t1"]



def test_extract_introns_weights_from_headers(tmp_path):
    ex = os.path.join(ROOT, "workflow", "scripts", "extract_introns.py")
    (tmp_path / "hq.fasta").write_text(">transcript/0 full_length_coverage=7;length=300;num_subreads=20\nACGT\n"
                                      ">transcript/1 full_length_coverage=2;length=300\nACGT\n>plain_read\nACGT\n")
    (tmp_path / "reads.fq").write_text("@transcript/1 full_length_coverage=4\nACGT\n+\nIIII\n")
    sam = ("transcript/0\t0\tc1\t100\t60\t50M200N50M\t*\t0\t0\tA\t*\n"
           "transcript/1\t0\tc1\t100\t60\t100M\t*\t0\t0\tA\t*\n"
           "plain_read\t0\tc1\t100\t60\t100M\t*\t0\t0\tA\t*\n")
    for reads, expect in (("hq.fasta", ["7", "2", "1"]), ("reads.fq", ["1", "4", "1"])):
        out = tmp_path / f"{reads}.tsv"
        subprocess.run([sys.executable, ex, "chains", "--reads", str(tmp_path / reads), "--out", str(out)],
                       input=sam, capture_output=True, text=True, check=True)
        assert [l.split("\t")[-1] for l in out.read_text().splitlines()[1:]] == expect


def test_r3_single_hq_transcript_with_coverage_supports(tmp_path):
    # one HQ transcript built from 2 full-length reads -> weight 2 >= min_lr_support 2 -> replace
    reads = [read_line(4950, 5950, [(5301, 5600)], weight=2)]
    dec, _ = run(tmp_path, [gene("e3", EVM3, src="EVM")], [gene("v3", EVI3)], [(EVI3, "+")], reads)
    assert dec[0]["eviann_support"] == "2" and dec[0]["action"] == "replace_with_eviann"
    assert "weighted by full_length_coverage" in (tmp_path / "dec.tsv").read_text().splitlines()[0]
    # the same transcript without a weight column counts as 1 read -> not enough
    dec, _ = run(tmp_path, [gene("e3", EVM3, src="EVM")], [gene("v3", EVI3)], [(EVI3, "+")],
                 [read_line(4950, 5950, [(5301, 5600)])])
    assert dec[0]["eviann_support"] == "1" and dec[0]["action"] == "keep_evm"


def test_concatenated_sample_tables(tmp_path):
    # per-sample tables concatenated by lr_read_chains keep only the first header; a repeated header
    # line must also be tolerated
    reads = ["seqid\tstart\tend\tintrons\tweight\n", read_line(4950, 5950, [(5301, 5600)], weight=1),
             read_line(4950, 5950, [(5301, 5600)], weight=1)]
    dec, _ = run(tmp_path, [gene("e3", EVM3, src="EVM")], [gene("v3", EVI3)], [(EVI3, "+")], reads)
    assert dec[0]["eviann_support"] == "2"


# ------------------------------------------------------------------ short-read support (A, B, C)
def short_support(tmp_path, junctions, stringtie, **kw):
    sys.path.insert(0, os.path.join(ROOT, "workflow", "scripts"))
    from refine import ShortReadSupport
    (tmp_path / "j.tsv").write_text("seqid\tstart\tend\treads\n" +
                                    "".join(f"c1\t{s}\t{e}\t{n}\n" for s, e, n in junctions))
    lines = []
    for i, (strand, cov, exons) in enumerate(stringtie):
        lines.append(f'c1\tStringTie\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t'
                     f'transcript_id "S{i}"; cov "{cov}";\n')
        lines += [f'c1\tStringTie\texon\t{s}\t{e}\t.\t{strand}\t.\ttranscript_id "S{i}";\n' for s, e in exons]
    (tmp_path / "st.gtf").write_text("".join(lines))
    return ShortReadSupport(str(tmp_path / "j.tsv"), str(tmp_path / "st.gtf"), kw.get("min_reads", 3),
                            kw.get("ratio", 0.1), kw.get("cov", 3), 0.5)


def test_short_support_b_absolute_and_relative(tmp_path):
    ex = [(100, 200), (301, 400)]                                   # intron 201-300
    assert short_support(tmp_path, [(201, 300, 10)], []).basis("c1", ex, "+") == "B"
    assert short_support(tmp_path, [(201, 300, 2)], []).basis("c1", ex, "+") == "-"        # < 3 reads
    # a competing overlapping intron with 200 reads: 10/200 = 0.05 < 0.1 -> not supported
    sup = short_support(tmp_path, [(201, 300, 10), (201, 350, 200)], [])
    assert sup.basis("c1", ex, "+") == "-" and sup.support("c1", ex, "+") == 0
    assert short_support(tmp_path, [(201, 300, 30), (201, 350, 200)], []).support("c1", ex, "+") == 30


def test_short_support_a_needs_stringtie_chain(tmp_path):
    ex = [(100, 200), (301, 400), (501, 600)]                       # introns 201-300, 401-500
    j = [(201, 300, 10), (401, 500, 12)]
    assert short_support(tmp_path, j, []).basis("c1", ex, "+") == "-"                      # B only
    chain = [(50, 200), (301, 400), (501, 700)]                     # same introns, longer UTRs: ok
    assert short_support(tmp_path, j, [("+", 5, chain)]).basis("c1", ex, "+") == "A+B"
    assert short_support(tmp_path, j, [("+", 2, chain)]).basis("c1", ex, "+") == "-"       # cov < 3
    assert short_support(tmp_path, j, [("-", 5, chain)]).basis("c1", ex, "+") == "-"       # other strand
    extra = [(100, 200), (301, 340), (361, 400), (501, 600)]       # extra intron inside the model
    assert short_support(tmp_path, j, [("+", 5, extra)]).basis("c1", ex, "+") == "-"
    assert short_support(tmp_path, j, [(".", 5, chain)]).support("c1", ex, "+") == 10      # weakest intron


def test_short_support_c_single_exon(tmp_path):
    ex = [(1000, 1600)]
    st = [("+", 4, [(950, 1500)])]                                  # covers 501/601 bp >= 50%
    assert short_support(tmp_path, [], st).basis("c1", ex, "+") == "C"
    assert short_support(tmp_path, [], [("+", 4, [(1400, 1700)])]).basis("c1", ex, "+") == "-"   # < 50%
    assert short_support(tmp_path, [(1100, 1300, 8)], st).basis("c1", ex, "+") == "-"            # spliced inside
    assert short_support(tmp_path, [(1100, 1300, 2)], st).basis("c1", ex, "+") == "C"            # weak junction


def test_no_rna_mode_all_r5(tmp_path):
    """Annotation without RNA-seq: empty EviAnn placeholder, no read evidence, empty --stringtie."""
    p = lambda n: str(tmp_path / n)
    open(p("g.fa"), "w").write(">c1\n" + "C" * L + "\n")
    open(p("evm.gff3"), "w").write(gene("prot", [(11001, 11300)], src="EVM") + gene("bare", [(12001, 12300)], src="EVM"))
    open(p("evi.gff3"), "w").write("##gff-version 3\n")
    open(p("mp.gff3"), "w").write("c1\tminiprot\tCDS\t11050\t11200\t.\t+\t0\tParent=MP1\n")
    open(p("te.gff3"), "w").write("")
    r = subprocess.run([sys.executable, SCRIPT, "--evm", p("evm.gff3"), "--eviann", p("evi.gff3"), "--genome", p("g.fa"),
                        "--stringtie", "", "--proteins", p("mp.gff3"), "--te-gff", p("te.gff3"),
                        "--out", p("out.gff3"), "--decisions", p("dec.tsv")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    lines = open(p("dec.tsv")).read().splitlines()
    assert lines[0].startswith("# support=none")
    dec = list(csv.DictReader(lines[1:], delimiter="\t"))
    assert {x["rule"] for x in dec} == {"R5"}
    assert {x["evm_genes"]: x["action"] for x in dec} == {"prot": "keep_evm", "bare": "keep_evm_flagged"}


def test_extract_introns_weights_from_count_tables(tmp_path):
    """FASTQ HQ transcripts without full_length_coverage: weights from an Iso-Seq count table."""
    ex = os.path.join(ROOT, "workflow", "scripts", "extract_introns.py")
    (tmp_path / "hq.fq").write_text("@transcript/0\nACGT\n+\nIIII\n@transcript/1\nACGT\n+\nIIII\n")   # no tags
    sam = ("transcript/0\t0\tc1\t100\t60\t50M200N50M\t*\t0\t0\tA\t*\n"
           "transcript/1\t0\tc1\t100\t60\t100M\t*\t0\t0\tA\t*\n"
           "transcript/2\t0\tc1\t100\t60\t100M\t*\t0\t0\tA\t*\n")
    tables = {
        "cluster_report.csv": "cluster_id,read_id,read_type\ntranscript/0,m1/1/ccs,FL\ntranscript/0,m1/2/ccs,FL\n"
                              "transcript/0,m1/3/ccs,FL\ntranscript/1,m1/4/ccs,FL\n",
        "x.flnc_count.txt": "id,leaf,root\ntranscript/0,2,1\ntranscript/1,1,0\n",       # per-sample columns summed
        "simple.tsv": "transcript/0\t3\ntranscript/1\t1\n",                               # no header
    }
    for name, text in tables.items():
        (tmp_path / name).write_text(text)
        out, summ = tmp_path / f"{name}.chains", tmp_path / f"{name}.summary"
        r = subprocess.run([sys.executable, ex, "chains", "--reads", str(tmp_path / "hq.fq"),
                            "--counts", str(tmp_path / name), "--sample", "iso1", "--platform", "pacbio_isoseq",
                            "--summary", str(summ), "--out", str(out)], input=sam, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        assert [l.split("\t")[-1] for l in out.read_text().splitlines()[1:]] == ["3", "1", "1"], name
        row = summ.read_text().splitlines()[1].split("\t")
        assert row[:6] == ["iso1", "pacbio_isoseq", "3", "2", "66.7", "5"] and name in row[6], name
    # headers without tags and no count table: everything weighs 1, summary shows 0 weighted reads
    summ = tmp_path / "none.summary"
    subprocess.run([sys.executable, ex, "chains", "--reads", str(tmp_path / "hq.fq"), "--sample", "iso1",
                    "--platform", "pacbio_isoseq", "--summary", str(summ), "--out", str(tmp_path / "n.tsv")],
                   input=sam, capture_output=True, text=True, check=True)
    assert summ.read_text().splitlines()[1].split("\t")[3] == "0"
