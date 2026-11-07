"""rename_genes.py: chromosome-based gene names."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "rename_genes.py")

# internal names seq0001.. were produced by the renaming step; originals are in the map
MAP = "new_id\toriginal_id\nseq0001\tChr1\nseq0002\tChr2\nseq0003\tscaf15\nseq0004\tscaf16\nseq0005\tChrM\n"


def gene(seq, gid, strand, txs, src="EVM", extra=""):
    lo = min(e[0] for _, ex, _ in txs for e in ex)
    hi = max(e[1] for _, ex, _ in txs for e in ex)
    out = [f"{seq}\t{src}\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID={gid}{extra}"]
    for tid, ex, cds in txs:
        out.append(f"{seq}\t{src}\tmRNA\t{ex[0][0]}\t{ex[-1][1]}\t.\t{strand}\t.\tID={tid};Parent={gid}")
        out += [f"{seq}\t{src}\texon\t{s}\t{e}\t.\t{strand}\t.\tParent={tid}" for s, e in ex]
        out += [f"{seq}\t{src}\tCDS\t{s}\t{e}\t.\t{strand}\t0\tParent={tid}" for s, e in cds]
    return "\n".join(out) + "\n"


CODING = (
    gene("seq0002", "b", "+", [("b.t", [(500, 900)], [(500, 799)])]) +
    gene("seq0001", "a2", "-", [("a2.t", [(5000, 5400), (5600, 6000)], [(5100, 5400), (5600, 5799)])],
         extra=";refine_rule=R1;eviann_support=yes") +
    gene("seq0001", "a1", "+", [("a1.short", [(100, 400)], [(150, 299)]),
                                ("a1.long", [(100, 400), (600, 1000)], [(150, 400), (600, 898)])]) +
    gene("seq0004", "s2", "+", [("s2.t", [(10, 400)], [(10, 309)])]) +
    gene("seq0003", "s1", "+", [("s1.t", [(10, 400)], [(10, 309)])]) +
    gene("seq0005", "m1", "+", [("m1.t", [(10, 400)], [(10, 309)])])
)
LNC = ("seq0001\tEviAnn\tgene\t3000\t3500\t.\t+\t.\tID=n1\n"
       "seq0001\tEviAnn\tlnc_RNA\t3000\t3500\t.\t+\t.\tID=n1.t;Parent=n1\n"
       "seq0001\tEviAnn\texon\t3000\t3500\t.\t+\t.\tParent=n1.t\n")


def run(tmp_path, *extra, coding=CODING, check=True):
    for n, t in (("c.gff3", coding), ("l.gff3", LNC), ("map.tsv", MAP)):
        (tmp_path / n).write_text(t)
    args = [sys.executable, SCRIPT, "--coding", tmp_path / "c.gff3", "--lncrna", tmp_path / "l.gff3",
            "--seqid-map", tmp_path / "map.tsv", "--prefix", "Mol", "--out-coding", tmp_path / "o.gff3",
            "--out-lncrna", tmp_path / "ol.gff3", "--out-map", tmp_path / "ids.tsv", *extra]
    r = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    if check:
        assert r.returncode == 0, r.stderr
    return r


def rows(path):
    return [l.split("\t") for l in path.read_text().splitlines() if not l.startswith("#")]


def ids_of(path, typ):
    return [dict(kv.split("=", 1) for kv in r[8].split(";"))["ID"] for r in rows(path) if r[2] == typ]


def test_chromosome_names_and_order(tmp_path):
    run(tmp_path)
    # genome order: Chr1 (a1, a2 by position), Chr2, scaf15, scaf16, ChrM (no number -> unplaced)
    assert ids_of(tmp_path / "o.gff3", "gene") == ["Mol01G000010", "Mol01G000020", "Mol02G000010",
                                                   "Mol00G000010", "Mol00G000020", "Mol00G000030"]
    m = {r["old_id"]: r for r in csv.DictReader(open(tmp_path / "ids.tsv"), delimiter="\t")}
    assert m["a1"]["new_id"] == "Mol01G000010" and m["a1"]["seqid"] == "Chr1"
    assert m["s1"]["new_id"] == "Mol00G000010" and m["s2"]["new_id"] == "Mol00G000020"
    assert ids_of(tmp_path / "ol.gff3", "gene") == ["Mol01L000010"]
    assert ids_of(tmp_path / "ol.gff3", "lnc_RNA") == ["Mol01L000010.1"]


def test_isoforms_attributes_utrs_phases(tmp_path):
    run(tmp_path)
    r = rows(tmp_path / "o.gff3")
    # .1 = longest CDS (a1.long)
    m = {x["old_id"]: x["new_id"] for x in csv.DictReader(open(tmp_path / "ids.tsv"), delimiter="\t")}
    assert m["a1.long"] == "Mol01G000010.1" and m["a1.short"] == "Mol01G000010.2"
    g2 = next(x for x in r if x[2] == "gene" and "ID=Mol01G000020;" in x[8])
    assert "refine_rule=R1" in g2[8] and "eviann_support=yes" in g2[8] and "source_id=a2" in g2[8]
    utr = {(x[2], x[3], x[4]) for x in r if "Mol01G000020.1" in x[8] and "UTR" in x[2]}
    # minus strand: genomic right of the CDS (5800-6000) is the 5' UTR, left (5000-5099) the 3' UTR
    assert utr == {("five_prime_UTR", "5800", "6000"), ("three_prime_UTR", "5000", "5099")}
    ph = {(x[3], x[4]): x[7] for x in r if x[2] == "CDS" and "Mol01G000020.1" in x[8]}
    assert ph[("5600", "5799")] == "0" and ph[("5100", "5400")] == "1"      # 200 bp first -> phase 1


def test_chrom_code_override_and_sequential(tmp_path):
    run(tmp_path, "--code", "ChrM=Mt")
    assert ids_of(tmp_path / "o.gff3", "gene")[-1] == "MolMtG000010"
    run(tmp_path, "--scheme", "sequential", "--step", "5", "--digits", "5")
    assert ids_of(tmp_path / "o.gff3", "gene")[:3] == ["MolG00005", "MolG00010", "MolG00015"]


def test_duplicate_chromosome_number_is_an_error(tmp_path):
    (tmp_path / "map.tsv").write_text("x")
    r = run(tmp_path, check=False, coding=CODING)
    assert r.returncode == 0                                      # sanity: default map is fine
    bad = MAP.replace("seq0002\tChr2", "seq0002\tchr01")           # Chr1 and chr01 -> both 1
    (tmp_path / "bad.tsv").write_text(bad)
    r = subprocess.run([sys.executable, SCRIPT, "--coding", str(tmp_path / "c.gff3"), "--seqid-map",
                        str(tmp_path / "bad.tsv"), "--prefix", "Mol", "--out-coding", str(tmp_path / "x.gff3"),
                        "--out-map", str(tmp_path / "x.tsv")], capture_output=True, text=True)
    assert r.returncode != 0 and "both map to chromosome 1" in r.stderr + r.stdout
