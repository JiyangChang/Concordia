"""repeat_stats.py on realistic EDTA / RepeatMasker / LTR_retriever output."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "repeat_stats.py")

TBL = """==================================================
file name: genome.fa
sequences:            2
total length:      10000 bp  (10000 bp excl N/X-runs)
GC level:         36.06 %
bases masked:       4500 bp ( 45.00 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Simple repeats:         3          150 bp    1.50 %
Low complexity:         2           50 bp    0.50 %
"""
# intact LTR: repeat_region parent + LTR_retrotransposon + 2 LTRs + TSD -> only the element counts
GFF = """##gff-version 3
c1\tEDTA\trepeat_region\t100\t1100\t.\t+\t.\tID=rr1;Classification=LTR/Gypsy;Method=structural
c1\tEDTA\ttarget_site_duplication\t100\t104\t.\t+\t.\tID=t1;Parent=rr1
c1\tEDTA\tlong_terminal_repeat\t105\t300\t.\t+\t.\tID=l1;Parent=rr1;Classification=LTR/Gypsy
c1\tEDTA\tLTR_retrotransposon\t105\t1095\t.\t+\t.\tID=e1;Parent=rr1;Classification=LTR/Gypsy;Method=structural
c1\tEDTA\tlong_terminal_repeat\t900\t1095\t.\t+\t.\tID=l2;Parent=rr1;Classification=LTR/Gypsy
c1\tEDTA\tCopia_LTR_retrotransposon\t2000\t2499\t.\t+\t.\tID=e2;Classification=LTR/Copia;Method=homology
c1\tEDTA\tCopia_LTR_retrotransposon\t2400\t2699\t.\t+\t.\tID=e3;Classification=LTR/Copia;Method=homology
c2\tEDTA\tTIR_transposon\t500\t699\t.\t-\t.\tID=e4;Classification=DNA/DTA;Method=homology
"""
LIB = ">TE_1#LTR/Gypsy\nACGT\n>TE_2#LTR/Gypsy\nACGT\n>TE_3#DNA/DTA\nACGT\n"
DIVSUM = """Weighted average Kimura divergence for each repeat family
Class\tRepeat\tabsLen\twellCharLen\tKimura%
LTR/Gypsy\tTE_1\t991\t991\t5.0

Coverage for each repeat class and divergence (Kimura)
Div LTR/Gypsy LTR/Copia DNA/DTA
0 100 50 0
1 300 0 200
"""
PASS = ("#LTR_loc\tCategory\tMotif\tTSD\t5'_TSD\t3'_TSD\tInternal\tIdentity\tStrand\tSuperFamily\tTE_type\tInsertion_Time\n"
        "c1:105..1095\tpass\tmotif:TGCA\tTSD:ACGTA\t100..104\t1096..1100\tIN:301..899\t0.9900\t+\tGypsy\tLTR\t384615\n")


def rows(path):
    return list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))


def test_repeat_stats(tmp_path):
    for n, t in (("te.gff3", GFF), ("g.tbl", TBL), ("lib.fa", LIB), ("g.divsum", DIVSUM), ("p.list", PASS)):
        (tmp_path / n).write_text(t)
    out = tmp_path / "stats"
    r = subprocess.run([sys.executable, SCRIPT, "--te-gff", tmp_path / "te.gff3", "--tbl", tmp_path / "g.tbl",
                        "--lib", tmp_path / "lib.fa", "--divsum", tmp_path / "g.divsum", "--pass-list", tmp_path / "p.list",
                        "--mutation-rate", "1.3e-8 (EDTA default, rice)", "--outdir", out], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    ov = {x["item"]: x for x in rows(out / "overview.tsv")}
    assert ov["genome size (bp)"]["bp"] == "10000"
    assert ov["masked by RepeatMasker (EDTA library, incl. simple/low complexity)"]["genome_%"] == "45.00"
    # Gypsy element 991 bp (children not double counted) + Copia 2000-2699 merged 700 bp + DTA 200 bp = 1891
    assert ov["TE annotated by EDTA (merged)"]["bp"] == "1891" and ov["TE annotated by EDTA (merged)"]["genome_%"] == "18.91"
    assert ov["Simple repeats (3 elements)"]["bp"] == "150"
    cl = {(x["class"], x["superfamily"]): x for x in rows(out / "te_classes.tsv")}
    assert cl[("LTR", "(all)")]["bp"] == "1691" and cl[("LTR", "(all)")]["copies"] == "3"
    assert cl[("LTR", "Gypsy")]["copies"] == "1" and cl[("LTR", "Gypsy")]["bp"] == "991"
    assert cl[("LTR", "Copia")]["copies"] == "2" and cl[("LTR", "Copia")]["bp"] == "700"
    assert cl[("DNA", "DTA")]["genome_%"] == "2.00" and cl[("total TE", "-")]["copies"] == "4"
    me = {x["method"]: x for x in rows(out / "te_methods.tsv")}
    assert me["structural"]["copies"] == "1" and me["homology"]["copies"] == "3"
    assert {x["classification"]: x["consensus_sequences"] for x in rows(out / "library.tsv")} == {"LTR/Gypsy": "2", "DNA/DTA": "1"}
    land = {(x["kimura_div"], x["class"]): x["bp"] for x in rows(out / "landscape.tsv")}
    assert land[("0", "LTR")] == "150" and land[("1", "LTR")] == "300" and land[("1", "DNA")] == "200"
    age = rows(out / "ltr_age.tsv")
    assert age[0]["superfamily"] == "Gypsy" and age[0]["insertion_mya"] == "0.385"
    assert "1.3e-8 (EDTA default, rice)" in (out / "ltr_age.tsv").read_text()


def test_repeat_stats_minimal_inputs(tmp_path):
    (tmp_path / "te.gff3").write_text(GFF)
    (tmp_path / "g.tbl").write_text(TBL)
    out = tmp_path / "stats"
    r = subprocess.run([sys.executable, SCRIPT, "--te-gff", tmp_path / "te.gff3", "--tbl", tmp_path / "g.tbl",
                        "--pass-list", "", "--outdir", out], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert rows(out / "landscape.tsv") == [] and rows(out / "ltr_age.tsv") == []


def test_pass_list_space_separated_header(tmp_path):
    # LTR_retriever writes the pass.list header with spaces
    space = ("#LTR_loc        Category        Motif   TSD     5'_TSD  3'_TSD  Internal        Identity        "
             "Strand  SuperFamily     TE_type Insertion_Time\n"
             "c1:105..1095\tpass\tmotif:TGCA\tTSD:ACGTA\t100..104\t1096..1100\tIN:301..899\t0.9900\t+\tCopia\tLTR\t769230\n")
    for n, t in (("te.gff3", GFF), ("g.tbl", TBL), ("p.list", space)):
        (tmp_path / n).write_text(t)
    out = tmp_path / "stats"
    r = subprocess.run([sys.executable, SCRIPT, "--te-gff", tmp_path / "te.gff3", "--tbl", tmp_path / "g.tbl",
                        "--pass-list", tmp_path / "p.list", "--outdir", out], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    age = rows(out / "ltr_age.tsv")
    assert age[0]["superfamily"] == "Copia" and age[0]["insertion_mya"] == "0.769"


def test_lowercase_edta_attributes(tmp_path):
    # EDTA writes classification= / method= in lowercase
    low = GFF.replace("Classification=", "classification=").replace("Method=", "method=")
    (tmp_path / "te.gff3").write_text(low)
    (tmp_path / "g.tbl").write_text(TBL)
    out = tmp_path / "stats"
    subprocess.run([sys.executable, SCRIPT, "--te-gff", tmp_path / "te.gff3", "--tbl", tmp_path / "g.tbl",
                    "--outdir", out], check=True, capture_output=True)
    cl = {(x["class"], x["superfamily"]): x["bp"] for x in rows(out / "te_classes.tsv")}
    assert cl[("LTR", "Gypsy")] == "991"
    assert {x["method"] for x in rows(out / "te_methods.tsv")} == {"structural", "homology"}
