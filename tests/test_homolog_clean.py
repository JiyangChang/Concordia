"""homolog_clean.py: keep protein-coding genes only."""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "homolog_clean.py")

# Araport-style typed genes + Ensembl/NCBI-style biotype attributes
GFF = """##gff-version 3
c1\ta\tprotein_coding_gene\t100\t900\t.\t+\t.\tID=PC1;Name=PC1
c1\ta\tmRNA\t100\t900\t.\t+\t.\tID=PC1.1;Parent=PC1
c1\ta\texon\t100\t900\t.\t+\t.\tParent=PC1.1
c1\ta\tfive_prime_UTR\t100\t199\t.\t+\t.\tParent=PC1.1
c1\ta\tCDS\t200\t799\t.\t+\t0\tID=PC1.1:CDS:1;Parent=PC1.1
c1\ta\tmRNA\t100\t900\t.\t+\t.\tID=PC1.2;Parent=PC1
c1\ta\texon\t100\t900\t.\t+\t.\tParent=PC1.2
c1\ta\ttransposable_element_gene\t1000\t2000\t.\t+\t.\tID=TE1
c1\ta\ttranscript\t1000\t2000\t.\t+\t.\tID=TE1.1;Parent=TE1
c1\ta\tCDS\t1000\t1299\t.\t+\t0\tParent=TE1.1
c1\ta\tpseudogene\t3000\t3500\t.\t+\t.\tID=PS1
c1\ta\tpseudogenic_transcript\t3000\t3500\t.\t+\t.\tID=PS1.1;Parent=PS1
c1\ta\tncRNA_gene\t4000\t4500\t.\t+\t.\tID=NC1
c1\ta\tlnc_RNA\t4000\t4500\t.\t+\t.\tID=NC1.1;Parent=NC1
c1\ta\texon\t4000\t4500\t.\t+\t.\tParent=NC1.1
c1\tb\tgene\t5000\t6000\t.\t-\t.\tID=G2;gene_biotype=protein_coding
c1\tb\tmRNA\t5000\t6000\t.\t-\t.\tID=G2.t;Parent=G2
c1\tb\tCDS\t5000\t5299\t.\t-\t0\tParent=G2.t
c1\tb\tgene\t7000\t8000\t.\t+\t.\tID=G3;gene_biotype=transposable_element
c1\tb\tmRNA\t7000\t8000\t.\t+\t.\tID=G3.t;Parent=G3
c1\tb\tCDS\t7000\t7299\t.\t+\t0\tParent=G3.t
c1\tb\tgene\t9000\t9900\t.\t+\t.\tID=G4;pseudo=true
c1\tb\tmRNA\t9000\t9900\t.\t+\t.\tID=G4.t;Parent=G4
c1\tb\tCDS\t9000\t9299\t.\t+\t0\tParent=G4.t
"""


def run(*args, check=True):
    r = subprocess.run([sys.executable, SCRIPT, *map(str, args)], capture_output=True, text=True)
    if check:
        assert r.returncode == 0, r.stderr
    return r


def clean(tmp_path, *extra):
    (tmp_path / "in.gff3").write_text(GFF)
    run("gff", "--gff", tmp_path / "in.gff3", "--out-gff", tmp_path / "out.gff3", "--out-ids", tmp_path / "ids",
        "--out-report", tmp_path / "rep.tsv", *extra)
    rows = [l.split("\t") for l in (tmp_path / "out.gff3").read_text().splitlines() if not l.startswith("#")]
    rep = dict(l.split("\t") for l in (tmp_path / "rep.tsv").read_text().splitlines()[1:])
    return rows, (tmp_path / "ids").read_text().split(), rep


def test_keeps_only_protein_coding(tmp_path):
    rows, ids, rep = clean(tmp_path)
    genes = [r[8].split(";")[0] for r in rows if r[2] == "gene"]
    assert genes == ["ID=PC1", "ID=G2;gene_biotype=protein_coding".split(";")[0]]
    assert ids == ["PC1.1", "G2.t"]                            # PC1.2 has no CDS -> dropped
    assert {r[2] for r in rows} == {"gene", "mRNA", "exon", "five_prime_UTR", "CDS"}
    assert sum(r[2] == "CDS" for r in rows) == 2
    assert rep["kept genes"] == "2" and rep["dropped gene: type transposable_element_gene"] == "1"
    assert rep["dropped gene: type pseudogene"] == "1" and rep["dropped gene: type ncRNA_gene"] == "1"
    assert rep["dropped gene: non-coding biotype / pseudo"] == "2"   # G3 (TE biotype), G4 (pseudo=true)
    assert rep["dropped transcripts of kept genes (no CDS / non-coding)"] == "1"


def test_passthrough(tmp_path):
    rows, ids, rep = clean(tmp_path, "--passthrough")
    assert (tmp_path / "out.gff3").read_text() == GFF
    assert "TE1.1" in ids and "PC1.2" not in ids
    assert any("disabled" in k for k in rep)


def test_fasta_selection_and_fallback(tmp_path):
    (tmp_path / "ids").write_text("PC1.1\nG2.t\n")
    (tmp_path / "p.fa").write_text(">PC1.1 | PC1\nMAAA*\n>G2.t\nMVV\n>TE1.1\nMKK\n")
    run("fasta", "--fasta", tmp_path / "p.fa", "--ids", tmp_path / "ids", "--out", tmp_path / "o.fa")
    assert (tmp_path / "o.fa").read_text() == ">PC1.1\nMAAA\n>G2.t\nMVV\n"
    (tmp_path / "q.fa").write_text(">other1\nMA\n>other2\nMV\n>PC1.1\nMA\n")
    r = run("fasta", "--fasta", tmp_path / "q.fa", "--ids", tmp_path / "ids", "--out", tmp_path / "o2.fa", check=False)
    assert r.returncode == 3                                   # 1/3 < 50% -> caller falls back to gffread
