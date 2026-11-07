"""Helper scripts in workflow/scripts/ (run as subprocesses, like the workflow does)."""
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "workflow", "scripts")


def run(name, *args):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, name), *map(str, args)],
                          capture_output=True, text=True)


def rename(tmp_path, fasta_text, **kw):
    fa = tmp_path / "g.fa"
    fa.write_text(fasta_text)
    out = [tmp_path / "o.fa", tmp_path / "map.tsv", tmp_path / "rep.tsv"]
    args = ["--fasta", fa, "--out-fasta", out[0], "--out-map", out[1], "--out-report", out[2]]
    for k, v in kw.items():
        args += ["--" + k.replace("_", "-"), v]
    r = run("rename_seqids.py", *args)
    return r, out


def test_rename_only_when_too_long(tmp_path):
    r, (fa, mp, rep) = rename(tmp_path, ">chr1\nACGT\n>chr2 desc\nAC\n")
    assert r.returncode == 0, r.stderr
    assert "renamed\tno" in rep.read_text()
    assert fa.read_text().startswith(">chr1\n")
    r, (fa, mp, rep) = rename(tmp_path, ">chr1\nACGT\n>scaffold_0000000001\nAC\n")
    assert "renamed\tyes" in rep.read_text()
    assert [l.split("\t")[0] for l in mp.read_text().splitlines()[1:]] == ["seq0001", "seq0002"]  # all renamed


def test_rename_link_unchanged(tmp_path):
    fa_in = tmp_path / "g.fa"
    # unchanged input -> symbolic link to the input, map still written
    r = run("rename_seqids.py", "--fasta", fa_in, "--out-fasta", tmp_path / "o.fa", "--out-map", tmp_path / "m.tsv",
            "--out-report", tmp_path / "r.tsv", "--link-unchanged") if fa_in.write_text(">chr1\nACGT\n") else None
    assert r.returncode == 0, r.stderr
    assert os.path.islink(tmp_path / "o.fa") and os.path.realpath(tmp_path / "o.fa") == str(fa_in.resolve())
    assert "symbolic link" in (tmp_path / "r.tsv").read_text()
    # renamed or filtered -> written file (an existing link is replaced, the input is not touched)
    fa_in.write_text(">scaffold_0000000001\nACGT\n")
    r = run("rename_seqids.py", "--fasta", fa_in, "--out-fasta", tmp_path / "o.fa", "--out-map", tmp_path / "m.tsv",
            "--out-report", tmp_path / "r.tsv", "--link-unchanged")
    assert r.returncode == 0, r.stderr
    assert not os.path.islink(tmp_path / "o.fa") and (tmp_path / "o.fa").read_text().startswith(">seq0001")
    assert fa_in.read_text().startswith(">scaffold_0000000001")


def test_rename_invalid_chars(tmp_path):
    r, (_, _, rep) = rename(tmp_path, ">chr1.1\nACGT\n")
    assert "renamed\tyes" in rep.read_text()
    r, _ = rename(tmp_path, ">chr1.1\nACGT\n", on_invalid_chars="false")
    assert r.returncode != 0 and "invalid characters" in r.stderr


def test_rename_all_filtered_is_error(tmp_path):
    r, _ = rename(tmp_path, ">chr1\nACGT\n", min_len=1000)
    assert r.returncode != 0 and "no sequences left" in r.stderr


def test_restore_roundtrip(tmp_path):
    mp = tmp_path / "map.tsv"
    mp.write_text("new_id\toriginal_id\nseq0001\tscaffold_0000000001\n")
    gff = tmp_path / "a.gff3"
    gff.write_text("##gff-version 3\n##sequence-region seq0001 1 100\nseq0001\tx\tgene\t1\t10\t.\t+\t.\tID=g1\n")
    out = tmp_path / "b.gff3"
    assert run("restore_seqids.py", "--gff", gff, "--map", mp, "--out", out).returncode == 0
    t = out.read_text()
    assert "scaffold_0000000001\tx\tgene" in t and "##sequence-region scaffold_0000000001 1 100" in t


GEMOMA_LIKE = """##gff-version 3
c1\tGeMoMa\tgene\t100\t400\t.\t+\t.\tID=g1
c1\tGeMoMa\tmRNA\t100\t400\t.\t+\t.\tID=t1;Parent=g1
c1\tGeMoMa\tCDS\t100\t200\t.\t+\t0\tParent=t1
c1\tGeMoMa\tCDS\t300\t400\t.\t+\t1\tParent=t1
c1\tx\tgene\t500\t900\t.\t+\t.\tID=nc1
c1\tx\tlnc_RNA\t500\t900\t.\t+\t.\tID=nct1;Parent=nc1
c1\tx\texon\t500\t900\t.\t+\t.\tParent=nct1
"""


def test_genepred_derives_exons_and_prefixes_ids(tmp_path):
    gff = tmp_path / "g.gff3"; gff.write_text(GEMOMA_LIKE)
    out = tmp_path / "o.gff3"
    r = run("evm_prepare.py", "genepred", "--gff", gff, "--source", "GeMoMa", "--out", out)
    assert r.returncode == 0, r.stderr
    rows = [l.split("\t") for l in out.read_text().splitlines()]
    assert [x[2] for x in rows] == ["gene", "mRNA", "exon", "exon", "CDS", "CDS"]  # non-coding dropped
    assert all(x[1] == "GeMoMa" for x in rows)
    assert rows[0][8] == "ID=GeMoMa.g1" and rows[2][8].endswith("Parent=GeMoMa.t1")


def test_split_coding(tmp_path):
    gff = tmp_path / "g.gff3"; gff.write_text(GEMOMA_LIKE)
    c, n = tmp_path / "c.gff3", tmp_path / "n.gff3"
    assert run("split_coding.py", "--gff", gff, "--coding", c, "--noncoding", n).returncode == 0
    assert "ID=g1" in c.read_text() and "nc1" not in c.read_text()
    assert "ID=nc1" in n.read_text()


def test_transcripts_target_coords_minus_strand(tmp_path):
    gtf = tmp_path / "s.gtf"
    gtf.write_text('c1\tStringTie\ttranscript\t100\t400\t.\t-\t.\ttranscript_id "T1"; gene_id "G1";\n'
                   'c1\tStringTie\texon\t100\t199\t.\t-\t.\ttranscript_id "T1"; gene_id "G1";\n'
                   'c1\tStringTie\texon\t301\t400\t.\t-\t.\ttranscript_id "T1"; gene_id "G1";\n')
    out = tmp_path / "t.gff3"
    assert run("evm_prepare.py", "transcripts", "--gtf", gtf, "--source", "LR_stringtie", "--out", out).returncode == 0
    rows = [l.split("\t") for l in out.read_text().splitlines()]
    # minus strand: the 3'-most genomic exon is transcript position 1
    assert rows[0][3] == "301" and rows[0][8].endswith("Target=T1 1 100 +")
    assert rows[1][3] == "100" and rows[1][8].endswith("Target=T1 101 200 +")


def test_miniprot_conversion(tmp_path):
    g = tmp_path / "m.gff3"
    g.write_text("c1\tminiprot\tmRNA\t100\t400\t50\t+\t.\tID=MP1;Target=P1 1 60\n"
                 "c1\tminiprot\tCDS\t100\t190\t50\t+\t0\tParent=MP1;Target=P1 1 30\n")
    out = tmp_path / "p.gff3"
    assert run("evm_prepare.py", "miniprot", "--gff", g, "--out", out).returncode == 0
    row = out.read_text().split("\t")
    assert row[2] == "nucleotide_to_protein_match" and "ID=miniprot.MP1;Target=P1 1 30" in row[8]


def test_qc_summary_parses_busco(tmp_path):
    d = tmp_path / "qc"
    (d / "busco" / "lin").mkdir(parents=True)
    (d / "busco" / "lin" / "short_summary.txt").write_text("\tC:95.1%[S:90.0%,D:5.1%],F:1.0%,M:3.9%,n:1614\n")
    out = tmp_path / "s.tsv"
    assert run("qc_summary.py", "--out", out, "--lineages", "lin", "--", f"final={d}").returncode == 0
    row = out.read_text().splitlines()[1].split("\t")
    assert row[:8] == ["final", "lin", "95.1", "90.0", "5.1", "1.0", "3.9", "1614"]


def test_split_coding_eviann_pseudogenes(tmp_path):
    gff = tmp_path / "e.gff3"
    gff.write_text(GEMOMA_LIKE +
                   "c1\tEviAnn\tgene\t9000\t9900\t.\t+\t.\tID=ps1;gene_biotype=processed_pseudogene;pseudo=true\n"
                   "c1\tEviAnn\tmRNA\t9000\t9900\t.\t+\t.\tID=ps1.t;Parent=ps1\n"
                   "c1\tEviAnn\tCDS\t9000\t9299\t.\t+\t0\tParent=ps1.t\n")
    c, n, p = tmp_path / "c.gff3", tmp_path / "n.gff3", tmp_path / "p.gff3"
    r = run("split_coding.py", "--gff", gff, "--coding", c, "--noncoding", n, "--pseudogenes", p)
    assert r.returncode == 0, r.stderr
    assert "ps1" not in c.read_text() and "ID=ps1" in p.read_text() and "pseudogenes: 1" in r.stdout


def test_miniprot_stats(tmp_path):
    faa = tmp_path / "p.faa"
    faa.write_text(">P1\nMA\n>P2\nMC\n>P3\nMD\n>P4\nME\n")
    gff = tmp_path / "mp.gff3"
    gff.write_text("##gff-version 3\n"
                   "c1\tminiprot\tmRNA\t1\t90\t9\t+\t.\tID=ath_MP1;Identity=0.95;Target=P1 1 30\n"
                   "c1\tminiprot\tCDS\t1\t90\t9\t+\t0\tParent=ath_MP1;Target=P1 1 30\n"
                   "c1\tminiprot\tmRNA\t300\t390\t9\t+\t.\tID=ath_MP2;Identity=0.5;Frameshift=1;Target=P1 1 30\n"
                   "c1\tminiprot\tmRNA\t500\t590\t9\t-\t.\tID=ath_MP3;Identity=0.7;Target=P2 1 30\n")
    out = tmp_path / "s.tsv"
    r = run("miniprot_stats.py", "--species", f"ath={faa},{gff}", "--out", out)
    assert r.returncode == 0, r.stderr
    row = out.read_text().splitlines()[1].split("\t")
    assert row == ["ath", "4", "2", "50.0", "3", "0.717", "1", "1"]


def test_orthodb_subset(tmp_path):
    import gzip
    fa = tmp_path / "odb.fa.gz"
    with gzip.open(fa, "wt") as fh:
        fh.write(">10_0:1\tx\nMA\nMA\n>20_0:1\nMC\n>30_0:1\nMD\n>10_0:2\nME\n")
    lv = tmp_path / "levels.tab"
    lv.write_text("33090\tViridiplantae\t100\t10\t3\n3699\tBrassicales\t50\t5\t1\n41938\tMalvales\t20\t2\t1\n")
    l2s = tmp_path / "l2s.tab"
    l2s.write_text("2759\t10_0\t1\t{2759,33090,3699,10}\n2759\t20_0\t1\t{2759,33090,41938,20}\n"
                   "2759\t30_0\t1\t{2759,33090,30}\n")
    out, rep = tmp_path / "s.faa", tmp_path / "r.tsv"
    r = run("orthodb_subset.py", "--fasta", fa, "--levels", lv, "--level2species", l2s,
            "--select", "brassicales", "41938", "--out", out, "--report", rep)
    assert r.returncode == 0, r.stderr
    assert out.read_text() == ">10_0:1\tx\nMA\nMA\n>20_0:1\nMC\n>10_0:2\nME\n"   # multi-line records kept whole
    assert "proteins_in_subset\t3" in rep.read_text() and "species_in_subset\t2" in rep.read_text()
    r = run("orthodb_subset.py", "--fasta", fa, "--levels", lv, "--level2species", l2s, "--select", "Poales",
            "--out", out)
    assert r.returncode != 0 and "not found" in r.stderr
    r = run("orthodb_subset.py", "--fasta", fa, "--levels", lv, "--level2species", l2s, "--list")
    assert "3699\tBrassicales\t1\t2" not in r.stdout          # levels with < 2 species are not listed
    assert "33090\tViridiplantae\t3\t4" in r.stdout


def test_select_longest(tmp_path):
    g = tmp_path / "all.gff3"
    g.write_text("##gff-version 3\n"
                 "c1\tx\tgene\t1\t1000\t.\t+\t.\tID=g1\n"
                 "c1\tx\tmRNA\t100\t900\t.\t+\t.\tID=g1.1;Parent=g1\n"
                 "c1\tx\texon\t100\t900\t.\t+\t.\tParent=g1.1\n"
                 "c1\tx\tCDS\t200\t400\t.\t+\t0\tParent=g1.1\n"
                 "c1\tx\tmRNA\t1\t1000\t.\t+\t.\tID=g1.2;Parent=g1\n"     # same CDS, longer transcript
                 "c1\tx\texon\t1\t1000\t.\t+\t.\tParent=g1.2\n"
                 "c1\tx\tCDS\t200\t400\t.\t+\t0\tParent=g1.2\n"
                 "c1\tx\tmRNA\t150\t950\t.\t+\t.\tID=g1.3;Parent=g1\n"   # shorter CDS
                 "c1\tx\texon\t150\t950\t.\t+\t.\tParent=g1.3,g1.2\n"
                 "c1\tx\tCDS\t300\t400\t.\t+\t0\tParent=g1.3\n"
                 "c1\tx\tgene\t2000\t3000\t.\t-\t.\tID=g2\n"
                 "c1\tx\tmRNA\t2100\t2900\t.\t-\t.\tID=g2.1;Parent=g2\n"
                 "c1\tx\tCDS\t2200\t2800\t.\t-\t0\tParent=g2.1\n")
    out, rep = tmp_path / "l.gff3", tmp_path / "r.tsv"
    r = run("select_longest.py", "--gff", g, "--out", out, "--report", rep)
    assert r.returncode == 0, r.stderr
    rows = [l.split("\t") for l in out.read_text().splitlines() if not l.startswith("#")]
    assert [x[8] for x in rows if x[2] == "mRNA"] == ["ID=g1.2;Parent=g1", "ID=g2.1;Parent=g2"]
    assert [x for x in rows if x[8] == "ID=g1"][0][3:5] == ["1", "1000"]          # gene = kept transcript
    assert not any("g1.1" in x[8] or "g1.3" in x[8] for x in rows)
    assert ["150", "950", "Parent=g1.2"] in [[x[3], x[4], x[8]] for x in rows]      # shared exon: parent trimmed
    assert "transcripts_kept\t2" in rep.read_text() and "genes_with_several_isoforms\t1" in rep.read_text()


def test_internal_stops_marked_and_cleaned(tmp_path):
    pep = tmp_path / "p.fa"
    pep.write_text(">t1\nMAK.\n>t2\nMA.K*\n>t3\n*\n")
    out = tmp_path / "c.fa"
    r = run("clean_pep.py", pep, out)
    assert r.returncode == 0 and out.read_text() == ">t1\nMAK\n>t2\nMAXK\n" and "1 protein" in r.stderr
    gff = tmp_path / "a.gff3"
    gff.write_text("##gff-version 3\nc1\tx\tgene\t1\t9\t.\t+\t.\tID=g2\n"
                   "c1\tx\tmRNA\t1\t9\t.\t+\t.\tID=t2;Parent=g2\nc1\tx\tCDS\t1\t9\t.\t+\t0\tParent=t2\n"
                   "c1\tx\tgene\t20\t29\t.\t+\t.\tID=g1\nc1\tx\tmRNA\t20\t29\t.\t+\t.\tID=t1;Parent=g1\n")
    rep = tmp_path / "s.tsv"
    r = run("mark_internal_stops.py", "--gff", gff, "--pep", pep, "--report", rep)
    assert r.returncode == 0, r.stderr
    rows = [l.split("\t") for l in gff.read_text().splitlines() if not l.startswith("#")]
    marked = [x[8] for x in rows if "internal_stop=yes" in x[8]]
    assert marked == ["ID=g2;internal_stop=yes", "ID=t2;Parent=g2;internal_stop=yes"]
    assert rep.read_text().splitlines()[1].startswith("t2\t")


def test_merge_function_any_excludes_name_fields(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text("".join(f"c1\tx\tgene\t{i}\t{i + 9}\t.\t+\t.\tID=g{i}\nc1\tx\tmRNA\t{i}\t{i + 9}\t.\t+\t.\tID=t{i};Parent=g{i}\n"
                           for i in (1, 2, 3)))
    cols = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
            "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway"]
    rows = [["t1", "x", "0", "1", "-", "-", "S", "Some protein", "ABC1", "-", "-", "-", "-"],   # COG -> any
            ["t2", "x", "0", "1", "-", "-", "-", "Only a description", "-", "-", "-", "-", "-"]]  # name only
    egg = tmp_path / "e.annotations"
    egg.write_text("## eggNOG\n#" + "\t".join(cols) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    o = {k: tmp_path / f"{k}.out" for k in ("gff", "tsv", "go", "kegg", "stats")}
    r = run("merge_function.py", "--gff", gff, "--eggnog", egg, "--out-gff", o["gff"], "--out-tsv", o["tsv"],
            "--out-go", o["go"], "--out-kegg", o["kegg"], "--out-stats", o["stats"])
    assert r.returncode == 0, r.stderr
    st = {l.split("\t")[0]: l.split("\t")[1] for l in o["stats"].read_text().splitlines()[1:]}
    assert st["any"] == "1" and st["Description"] == "2" and st["Preferred_name"] == "1"
    assert "product=Only a description" in o["gff"].read_text()          # still written to the GFF3


def test_braker4_inputs_share_busco_downloads(tmp_path):
    tpl = tmp_path / "config.ini.example"
    tpl.write_text("[paths]\nsamples_file = samples.csv\n\n[PARAMS]\nskip_busco = 0\n")
    base = ["--sample", "AB", "--genome", "g.fa", "--masked", "m.fa", "--proteins", "p.fa",
            "--busco-lineage", "brassicales_odb12", "--template", tpl,
            "--out-csv", tmp_path / "s.csv", "--out-ini", tmp_path / "c.ini"]
    r = run("make_braker4_inputs.py", *base, "--busco-download-path", tmp_path / "db" / "busco_downloads")
    assert r.returncode == 0, r.stderr
    ini = (tmp_path / "c.ini").read_text()
    assert f"busco_download_path = {tmp_path / 'db' / 'busco_downloads'}" in ini
    assert f"compleasm_download_path = {tmp_path / 'db' / 'busco_downloads' / 'lineages'}" in ini
    assert f"samples_file = {tmp_path / 's.csv'}" in ini and "skip_busco = 1" in ini
    r = run("make_braker4_inputs.py", *base)                      # without the option: BRAKER4's default path
    assert r.returncode == 0 and "busco_download_path" not in (tmp_path / "c.ini").read_text()
