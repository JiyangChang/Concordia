"""annotation_stats.py: gene-model statistics used for the EviAnn report section."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "annotation_stats.py")

CODING = """##gff-version 3
c1\tEviAnn\tgene\t100\t1000\t.\t+\t.\tID=g1
c1\tEviAnn\tmRNA\t100\t1000\t.\t+\t.\tID=g1.t1;Parent=g1
c1\tEviAnn\texon\t100\t400\t.\t+\t.\tParent=g1.t1
c1\tEviAnn\texon\t601\t1000\t.\t+\t.\tParent=g1.t1
c1\tEviAnn\tCDS\t200\t400\t.\t+\t0\tParent=g1.t1
c1\tEviAnn\tCDS\t601\t900\t.\t+\t1\tParent=g1.t1
c1\tEviAnn\tmRNA\t100\t1000\t.\t+\t.\tID=g1.t2;Parent=g1
c1\tEviAnn\texon\t100\t1000\t.\t+\t.\tParent=g1.t2
c1\tEviAnn\tfive_prime_UTR\t100\t199\t.\t+\t.\tParent=g1.t2
c1\tEviAnn\tCDS\t200\t898\t.\t+\t0\tParent=g1.t2
c1\tEviAnn\tgene\t2000\t2600\t.\t-\t.\tID=g2
c1\tEviAnn\tmRNA\t2000\t2600\t.\t-\t.\tID=g2.t1;Parent=g2
c1\tEviAnn\texon\t2000\t2600\t.\t-\t.\tParent=g2.t1
c1\tEviAnn\tCDS\t2000\t2299\t.\t-\t0\tParent=g2.t1
"""
LNC = """##gff-version 3
c1\tEviAnn\tgene\t5000\t5800\t.\t+\t.\tID=n1
c1\tEviAnn\tlnc_RNA\t5000\t5800\t.\t+\t.\tID=n1.t1;Parent=n1
c1\tEviAnn\texon\t5000\t5800\t.\t+\t.\tParent=n1.t1
"""


def test_annotation_stats(tmp_path):
    (tmp_path / "c.gff3").write_text(CODING)
    (tmp_path / "n.gff3").write_text(LNC)
    out, dist = tmp_path / "s.tsv", tmp_path / "d.tsv"
    r = subprocess.run([sys.executable, SCRIPT, "--set", f"coding={tmp_path / 'c.gff3'}",
                        "--set", f"lncRNA={tmp_path / 'n.gff3'}", "--out", out, "--out-dist", dist],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    m = {row["metric"]: row for row in csv.DictReader(open(out), delimiter="\t")}
    c, n = (lambda k: m[k]["coding"]), (lambda k: m[k]["lncRNA"])
    assert c("genes") == "2" and c("transcripts") == "3" and c("transcripts_per_gene") == "1.50"
    assert c("multi_isoform_genes") == "1" and c("single_exon_transcripts") == "2"
    assert c("coding_transcripts") == "3"
    # g1.t1: derived 5'+3' UTR; g1.t2: explicit 5' UTR only (explicit features win); g2.t1 (minus): exon beyond
    # CDS on the genomic right = 5' UTR only
    assert c("coding_tx_with_5UTR_%") == "100.0"
    assert c("coding_tx_with_3UTR_%") == "33.3" and c("coding_tx_with_both_UTR_%") == "33.3"
    assert c("cds_length_aa_median") == "167"        # 167, 233, 100 aa
    assert n("genes") == "1" and n("coding_transcripts") == "0" and n("single_exon_transcripts_%") == "100.0"
    kinds = {(x["set"], x["kind"]) for x in csv.DictReader(open(dist), delimiter="\t")}
    assert ("coding", "cds_aa") in kinds and ("lncRNA", "isoforms_per_gene") in kinds
    assert ("lncRNA", "cds_aa") not in kinds
