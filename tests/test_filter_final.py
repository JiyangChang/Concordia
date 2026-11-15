"""filter_final.py: remove for TE reasons only if TE overlap > threshold AND no RNA evidence."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "filter_final.py")


def gene(gid, start, end, strand="+", source="EVM", cds=None, exons=None):
    """One gene/one transcript; cds/exons are lists of (start, end)."""
    cds = cds or [(start, end)]
    exons = exons or cds
    rows = [f"c1\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}",
            f"c1\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={gid}.t1;Parent={gid}"]
    rows += [f"c1\t{source}\texon\t{s}\t{e}\t.\t{strand}\t.\tParent={gid}.t1" for s, e in exons]
    rows += [f"c1\t{source}\tCDS\t{s}\t{e}\t.\t{strand}\t0\tParent={gid}.t1" for s, e in cds]
    return "\n".join(rows) + "\n"


def run(tmp_path, genes, te=(), stringtie=(), decisions=(), braker=(), support=None, tesorter=(), **kw):
    p = lambda n: str(tmp_path / n)
    open(p("in.gff3"), "w").write("".join(genes))
    open(p("te.gff3"), "w").write("".join(f"c1\tEDTA\tLTR_retrotransposon\t{s}\t{e}\t.\t+\t.\tID=te{i}\n"
                                          for i, (s, e) in enumerate(te)))
    open(p("st.gtf"), "w").write("".join(f'c1\tStringTie\texon\t{s}\t{e}\t.\t{sd}\t.\ttranscript_id "S{i}";\n'
                                         for i, (s, e, sd) in enumerate(stringtie)))
    with open(p("dec.tsv"), "w") as fh:
        fh.write("evm_gene\tevm_tx\teviann_gene\teviann_tx\tclass_code\tcds_aa\trule\taction\tnote\n")
        for g, code in decisions:
            fh.write(f"{g}\t{g}.t1\tv\tv.t\t{code}\t100\tR1\tkeep_evm\t\n")
    open(p("b4.gff3"), "w").write("".join(braker))
    open(p("b4.tsv"), "w").write(support or "")
    open(p("tes.tsv"), "w").write("#TE\tOrder\tSuperfamily\n" + "".join(f"{t}\tLTR\tCopia\n" for t in tesorter))
    args = [sys.executable, SCRIPT, "--gff", p("in.gff3"), "--te-gff", p("te.gff3"), "--transcripts", p("st.gtf"),
            "--refine-decisions", p("dec.tsv"), "--braker4-gff", p("b4.gff3"), "--braker4-support", p("b4.tsv"),
            "--tesorter", p("tes.tsv"), "--out", p("out.gff3"), "--removed", p("rm.tsv"), "--decisions", p("d.tsv")]
    for k, v in kw.items():
        args += ["--" + k.replace("_", "-")] + ([] if v is True else [str(v)])     # True = flag without value
    r = subprocess.run(args, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with open(p("d.tsv")) as fh:
        rows = {d["gene"]: d for d in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")}
    return rows, r.stdout


TE_GENE = gene("g1", 1001, 1600)          # 600 bp CDS
TE_80 = [(1001, 1480)]                    # 480/600 = 80 % of CDS


def test_te_overlap_without_rna_is_removed(tmp_path):
    d, _ = run(tmp_path, [TE_GENE], te=TE_80)
    assert d["g1"]["action"] == "remove" and d["g1"]["reason"] == "TE_overlap_0.80_no_RNA"


def test_te_overlap_below_threshold_is_kept(tmp_path):
    d, _ = run(tmp_path, [TE_GENE], te=[(1001, 1180)])  # 30 %
    assert d["g1"]["action"] == "keep"


def test_exactly_threshold_is_kept(tmp_path):
    d, _ = run(tmp_path, [TE_GENE], te=[(1001, 1300)])  # 50 %, rule is "> 50 %"
    assert d["g1"]["action"] == "keep"


def test_stringtie_support_rescues(tmp_path):
    d, _ = run(tmp_path, [TE_GENE], te=TE_80, stringtie=[(1500, 1700, "+")])
    assert d["g1"]["rna_stringtie"] == "yes" and d["g1"]["action"] == "keep"


def test_unstranded_stringtie_counts_opposite_strand_does_not(tmp_path):
    assert run(tmp_path, [TE_GENE], te=TE_80, stringtie=[(1500, 1700, ".")])[0]["g1"]["action"] == "keep"
    assert run(tmp_path, [TE_GENE], te=TE_80, stringtie=[(1500, 1700, "-")])[0]["g1"]["action"] == "remove"


def test_eviann_origin_rescues(tmp_path):
    d, _ = run(tmp_path, [gene("EviAnn.v1", 1001, 1600, source="EviAnn_refine")], te=TE_80)
    assert d["EviAnn.v1"]["rna_eviann"] == "yes" and d["EviAnn.v1"]["action"] == "keep"


def test_exact_eviann_match_rescues(tmp_path):
    assert run(tmp_path, [TE_GENE], te=TE_80, decisions=[("g1", "=")])[0]["g1"]["action"] == "keep"
    assert run(tmp_path, [TE_GENE], te=TE_80, decisions=[("g1", "j")])[0]["g1"]["action"] == "remove"


def test_braker4_support_rescues(tmp_path):
    b4 = gene("b1", 1100, 1500)
    d, out = run(tmp_path, [TE_GENE], te=TE_80, braker=[b4], support="gene_id\tRNA_support\nb1\t1\n")
    assert "RNA_support" in out
    assert d["g1"]["rna_braker4"] == "yes" and d["g1"]["action"] == "keep"
    d, _ = run(tmp_path, [TE_GENE], te=TE_80, braker=[b4], support="gene_id\tRNA_support\nb1\t0\n")
    assert d["g1"]["action"] == "remove"


def test_missing_gene_support_is_reported_not_fatal(tmp_path):
    d, out = run(tmp_path, [TE_GENE], te=TE_80, braker=[gene("b1", 1100, 1500)], support="")
    assert "criterion c not used" in out and d["g1"]["action"] == "remove"


def test_tesorter_hit_alone_does_not_remove(tmp_path):
    d, _ = run(tmp_path, [TE_GENE], tesorter=["g1.t1"])
    assert d["g1"]["tesorter"] == "LTR/Copia" and d["g1"]["action"] == "keep"


def test_intronic_te_does_not_count(tmp_path):
    g = gene("g1", 1001, 3000, cds=[(1001, 1300), (2701, 3000)])
    d, _ = run(tmp_path, [g], te=[(1301, 2700)])   # TE fills the whole intron
    assert d["g1"]["te_overlap_frac"] == "0.000" and d["g1"]["action"] == "keep"


def test_short_cds_is_removed_regardless_of_rna(tmp_path):
    d, _ = run(tmp_path, [gene("g1", 1001, 1090)], stringtie=[(1000, 1100, "+")])   # 30 aa
    assert d["g1"]["action"] == "remove" and d["g1"]["reason"].startswith("CDS_30aa")


def test_outputs_consistent(tmp_path):
    run(tmp_path, [TE_GENE, gene("g2", 5001, 5600)], te=TE_80)
    assert "ID=g2" in (tmp_path / "out.gff3").read_text() and "ID=g1" not in (tmp_path / "out.gff3").read_text()
    assert (tmp_path / "rm.tsv").read_text().splitlines()[1:] == ["g1\tTE_overlap_0.80_no_RNA"]


GS_HEADER = ("gene_id\ttranscript_id\tchrom\tstrand\tnum_introns\tintrons_sup_rnaseq\tpct_introns_sup_rnaseq\t"
             "introns_sup_protein\tpct_introns_sup_protein\tintrons_sup_any\tpct_introns_sup_any\tnum_exons\t"
             "exons_sup_rnaseq\tpct_exons_sup_rnaseq\texons_sup_protein\tpct_exons_sup_protein\texons_sup_any\t"
             "pct_exons_sup_any\n")


def test_braker4_gene_support_real_format(tmp_path):
    # BRAKER4 scripts/gene_support_summary.py: '#' summary lines, then header, one row per transcript
    support = ("# Gene support summary\n# transcripts: 3\n" + GS_HEADER +
               "b1\tb1.t1\tc1\t+\t2\t1\t50.0\t0\tNA\t1\t50.0\t3\t0\tNA\t0\tNA\t0\tNA\n"     # intron RNA support
               "b2\tb2.t1\tc1\t+\t2\t0\tNA\t2\t100.0\t2\t100.0\t3\t3\t100.0\t0\tNA\t3\t100.0\n"  # protein only
               "b3\tb3.t1\tc1\t+\t0\t0\tNA\t0\tNA\t0\tNA\t1\t1\t100.0\t0\tNA\t1\t100.0\n")     # single exon
    braker = gene("b1", 1100, 1500) + gene("b2", 3100, 3500) + gene("b3", 5100, 5500)
    genes = [TE_GENE, gene("g2", 3001, 3600), gene("g3", 5001, 5600)]
    te = [(1001, 1480), (3001, 3480), (5001, 5480)]
    d, out = run(tmp_path, genes, te=te, braker=[braker], support=support)
    assert "introns_sup_rnaseq" in out
    assert d["g1"]["rna_braker4"] == "yes" and d["g1"]["action"] == "keep"
    assert d["g2"]["rna_braker4"] == "no" and d["g2"]["action"] == "remove"        # protein support is not RNA
    assert d["g3"]["rna_braker4"] == "yes" and d["g3"]["action"] == "keep"        # single-exon: exon support


def test_protein_evidence_rescues_without_rna(tmp_path):
    (tmp_path / "mp.gff3").write_text("c1\tminiprot\tCDS\t1100\t1200\t.\t+\t0\tParent=MP1\n"
                                      "c1\tminiprot\tCDS\t5000\t5100\t.\t-\t0\tParent=MP2\n")
    other = gene("g2", 5001, 5600)                                   # protein hit only on the other strand
    d, _ = run(tmp_path, [TE_GENE, other], te=TE_80 + [(5001, 5480)], protein_evidence=tmp_path / "mp.gff3")
    assert d["g1"]["action"] == "keep" and d["g1"]["protein"] == "yes"
    assert d["g2"]["action"] == "remove" and d["g2"]["reason"] == "TE_overlap_0.80_no_RNA_no_protein"
    d, _ = run(tmp_path, [TE_GENE], te=TE_80)                        # without the option: unchanged
    assert d["g1"]["action"] == "remove" and "protein" not in d["g1"]


def test_short_genes_kept_with_evidence(tmp_path):
    short_rna = gene("s1", 1001, 1150)                 # 150 bp CDS = 50 aa
    short_none = gene("s2", 5001, 5150)
    genes, st = [short_rna, short_none], [(1001, 1150, "+")]
    d, _ = run(tmp_path, genes, stringtie=st, min_cds_aa=100)
    assert d["s1"]["action"] == "remove" and d["s2"]["action"] == "remove"          # default: length only
    d, _ = run(tmp_path, genes, stringtie=st, min_cds_aa=100, keep_short_supported=True)
    assert d["s1"]["action"] == "keep" and d["s1"]["reason"] == "short_CDS_50aa_supported"
    assert d["s2"]["action"] == "remove" and d["s2"]["reason"] == "CDS_50aa_lt_100"   # no RNA evidence
    out = open(tmp_path / "out.gff3").read()
    assert "ID=s1;short_cds=yes" in out and "s2" not in out
