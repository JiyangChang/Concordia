"""EVM weight tuning: region subsetting and grid scoring."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S = os.path.join(ROOT, "workflow", "scripts")


def run(name, *args):
    r = subprocess.run([sys.executable, os.path.join(S, name), *map(str, args)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


def tx(gid, cds, seq="s2", strand="+"):
    out = [f"{seq}\tEVM\tgene\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\tID={gid}",
           f"{seq}\tEVM\tmRNA\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\tID={gid}.t;Parent={gid}"]
    out += [f"{seq}\tEVM\texon\t{s}\t{e}\t.\t{strand}\t.\tParent={gid}.t" for s, e in cds]
    out += [f"{seq}\tEVM\tCDS\t{s}\t{e}\t.\t{strand}\t0\tParent={gid}.t" for s, e in cds]
    return "\n".join(out) + "\n"


def test_subset_auto_and_named(tmp_path):
    (tmp_path / "g.fa").write_text(">s1\nACGT\n>s2\nACGTACGTAC\n")
    (tmp_path / "map.tsv").write_text("new_id\toriginal_id\ns1\tChr1\ns2\tChr2\n")
    (tmp_path / "x.gff3").write_text("##gff-version 3\ns1\ta\tgene\t1\t4\t.\t+\t.\tID=a\ns2\ta\tgene\t1\t9\t.\t+\t.\tID=b\n")
    out = tmp_path / "o"
    run("subset_evm_inputs.py", "--genome", tmp_path / "g.fa", "--seqid-map", tmp_path / "map.tsv",
        "--outdir", out, tmp_path / "x.gff3")
    assert (out / "region.tsv").read_text().splitlines()[1] == "s2\tChr2\t10"          # longest
    assert (out / "x.gff3").read_text().splitlines()[1:] == ["s2\ta\tgene\t1\t9\t.\t+\t.\tID=b"]
    run("subset_evm_inputs.py", "--genome", tmp_path / "g.fa", "--seqid-map", tmp_path / "map.tsv",
        "--region", "Chr1", "--outdir", out, tmp_path / "x.gff3")
    assert (out / "region.tsv").read_text().splitlines()[1] == "s1\tChr1\t4"           # original name


def test_scoring_picks_evidence_consistent_combo(tmp_path):
    good = [(1001, 1300), (1501, 1800)]          # chain supported by reads and miniprot
    bad = [(1001, 1300), (1401, 1800)]           # different intron
    (tmp_path / "region.tsv").write_text("internal_id\toriginal_id\tlength\ns2\tChr2\t50000\n")
    (tmp_path / "lr.tsv").write_text("seqid\tstart\tend\tintrons\tweight\n"
                                     "s2\t950\t1850\t1301-1500\t3\n" "s2\t20000\t21000\t20101-20500\t2\n")
    (tmp_path / "mp.gff3").write_text("s2\tminiprot\tmRNA\t1001\t1800\t.\t+\t.\tID=MP1;Identity=0.95\n"
                                      "s2\tminiprot\tCDS\t1001\t1300\t.\t+\t0\tParent=MP1\n"
                                      "s2\tminiprot\tCDS\t1501\t1800\t.\t+\t0\tParent=MP1\n")
    for cid, model, extra in (("c00", bad, ""), ("c01", good, ""), ("c02", good, tx("g2", [(20001, 20100), (20501, 20600)]))):
        d = tmp_path / "grid" / cid
        d.mkdir(parents=True)
        (d / "weights.txt").write_text(f"ABINITIO_PREDICTION\tBRAKER4\t{cid[-1]}\n")
        (d / "tune.EVM.gff3").write_text(tx("g1", model) + extra)
    args = []
    for c in ("c00", "c01", "c02"):
        args += ["--combo", f"{c}={tmp_path / 'grid' / c / 'tune.EVM.gff3'}"]
    run("evm_tune_score.py", *args, "--weights-dir", tmp_path / "grid", "--region-file", tmp_path / "region.tsv",
        "--miniprot", tmp_path / "mp.gff3", "--lr-chains", tmp_path / "lr.tsv",
        "--out", tmp_path / "res.tsv", "--best", tmp_path / "best.txt")
    rows = list(csv.DictReader((l for l in open(tmp_path / "res.tsv") if not l.startswith("#")), delimiter="\t"))
    r = {x["combo"]: x for x in rows}
    assert r["c00"]["lr_chain"] == "0.0000" and r["c00"]["prot_chain"] == "0.0000" and r["c00"]["F"] == "0.0000"
    assert r["c01"]["P"] == "1.0000" and r["c01"]["lr_loci"] == "0.5000"            # misses the second locus
    # c02 also recovers the second read locus; its extra gene has RNA but no protein support and
    # still counts as confirmed (any evidence) -> P stays 1, R rises -> best
    assert r["c02"]["prot_chain"] == "0.5000" and r["c02"]["P"] == "1.0000"
    assert r["c02"]["lr_loci"] == "1.0000" and rows[0]["combo"] == "c02"
    assert (tmp_path / "best.txt").read_text() == "ABINITIO_PREDICTION\tBRAKER4\t2\n"
    assert "region Chr2 (s2, 50000 bp)" in (tmp_path / "res.tsv").read_text().splitlines()[0]


def test_scoring_short_reads_uses_stringtie_chain(tmp_path):
    three = [(1001, 1100), (1201, 1300), (1401, 1500)]          # introns 1101-1200, 1301-1400
    (tmp_path / "region.tsv").write_text("internal_id\toriginal_id\tlength\ns2\tChr2\t50000\n")
    (tmp_path / "j.tsv").write_text("seqid\tstart\tend\treads\ns2\t1101\t1200\t9\ns2\t1301\t1400\t9\n")
    (tmp_path / "mp.gff3").write_text("")
    d = tmp_path / "grid" / "c00"
    d.mkdir(parents=True)
    (d / "weights.txt").write_text("ABINITIO_PREDICTION\tBRAKER4\t1\n")
    (d / "tune.EVM.gff3").write_text(tx("g1", three))
    st = ('s2\tStringTie\ttranscript\t950\t1550\t.\t+\t.\ttranscript_id "S1"; cov "6";\n' +
          "".join(f's2\tStringTie\texon\t{s}\t{e}\t.\t+\t.\ttranscript_id "S1";\n'
                  for s, e in [(950, 1100), (1201, 1300), (1401, 1550)]))
    for name, gtf, expect in (("no_st", "", "0.0000"), ("st", st, "1.0000")):
        (tmp_path / f"{name}.gtf").write_text(gtf)
        run("evm_tune_score.py", "--combo", f"c00={d / 'tune.EVM.gff3'}", "--weights-dir", tmp_path / "grid",
            "--region-file", tmp_path / "region.tsv", "--miniprot", tmp_path / "mp.gff3",
            "--sr-junctions", tmp_path / "j.tsv", "--stringtie", tmp_path / f"{name}.gtf",
            "--out", tmp_path / f"{name}.tsv", "--best", tmp_path / "best.txt")
        row = list(csv.DictReader((l for l in open(tmp_path / f"{name}.tsv") if not l.startswith("#")),
                                  delimiter="\t"))[0]
        assert row["sr_chain"] == expect, name                  # 2 introns: junctions alone are not enough


def test_scoring_protein_only_without_rna(tmp_path):
    good, bad = [(1001, 1300), (1501, 1800)], [(1001, 1300), (1401, 1800)]
    (tmp_path / "region.tsv").write_text("internal_id\toriginal_id\tlength\ns2\tChr2\t50000\n")
    (tmp_path / "mp.gff3").write_text("s2\tminiprot\tmRNA\t1001\t1800\t.\t+\t.\tID=MP1;Identity=0.95\n"
                                      "s2\tminiprot\tCDS\t1001\t1300\t.\t+\t0\tParent=MP1\n"
                                      "s2\tminiprot\tCDS\t1501\t1800\t.\t+\t0\tParent=MP1\n")
    args = []
    for cid, model in (("c00", bad), ("c01", good)):
        d = tmp_path / "grid" / cid
        d.mkdir(parents=True)
        (d / "weights.txt").write_text(f"ABINITIO_PREDICTION\tBRAKER4\t{cid[-1]}\n")
        (d / "tune.EVM.gff3").write_text(tx("g1", model))
        args += ["--combo", f"{cid}={d / 'tune.EVM.gff3'}"]
    run("evm_tune_score.py", *args, "--weights-dir", tmp_path / "grid", "--region-file", tmp_path / "region.tsv",
        "--miniprot", tmp_path / "mp.gff3", "--stringtie", "", "--out", tmp_path / "res.tsv", "--best", tmp_path / "b.txt")
    lines = (tmp_path / "res.tsv").read_text().splitlines()
    assert "confirmed by any of prot_chain; R = mean(prot_loci)" in lines[0]
    rows = list(csv.DictReader(lines[1:], delimiter="\t"))
    assert rows[0]["combo"] == "c01" and rows[0]["P"] == "1.0000" and rows[0]["prot_loci"] == "1.0000"
