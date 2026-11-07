"""Genome assembly quality scripts and config toggles."""
import csv
import os
import subprocess
import sys

import pytest

from common import ConfigError
from test_common import make

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "workflow", "scripts")


def run(name, *args):
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, name), *map(str, args)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


def table(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_assembly_stats_levels_and_gaps(tmp_path):
    # s1: 40 + gap(10 N) + 50 ; s2: 30 with a 3-N run (not a gap at gap_min_n=10) ; s3: 20
    s1 = "A" * 20 + "G" * 20 + "N" * 10 + "C" * 50
    s2 = "AT" * 6 + "NNN" + "A" * 15
    fa = tmp_path / "g.fa"
    fa.write_text(f">s1\n{s1}\n>s2 desc\n{s2}\n>s3\n{'GC' * 10}\n")
    run("assembly_stats.py", "--fasta", fa, "--gap-min-n", 10, "--out-summary", tmp_path / "s.tsv", "--out-seqs", tmp_path / "q.tsv")
    m = {r["metric"]: r for r in table(tmp_path / "s.tsv")}
    assert m["count"]["scaffold"] == "3" and m["count"]["contig"] == "4"
    assert m["total_length"]["scaffold"] == "150" and m["total_length"]["contig"] == "140"
    assert m["N50"]["scaffold"] == "100" and m["L50"]["scaffold"] == "1"
    assert m["N50"]["contig"] == "40" and m["L50"]["contig"] == "2"   # contigs 50,40,30,20
    assert m["max_length"]["contig"] == "50"
    assert m["gaps"]["scaffold"] == "1" and m["gap_length"]["scaffold"] == "10"
    assert m["N_bases"]["scaffold"] == "13" and m["N_bases"]["contig"] == "3"
    seqs = table(tmp_path / "q.tsv")
    assert [s["seq_id"] for s in seqs] == ["s1", "s2", "s3"] and seqs[0]["contigs"] == "2"


def test_telomere_summary(tmp_path):
    (tmp_path / "seqs.tsv").write_text("seq_id\tlength\nchr1\t30000\nchr2\t30000\nchr3\t30000\n")
    (tmp_path / "w.tsv").write_text(
        "id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\n"
        "chr1\t10000\t0\t300\tTTTAGGG\nchr1\t20000\t0\t0\tTTTAGGG\nchr1\t30000\t250\t0\tTTTAGGG\n"
        "chr2\t10000\t0\t120\tTTTAGGG\nchr2\t30000\t5\t0\tTTTAGGG\n")
    run("telomere_summary.py", "--windows", tmp_path / "w.tsv", "--seqs", tmp_path / "seqs.tsv", "--motif", "TTTAGGG",
        "--min-repeats", 50, "--out", tmp_path / "t.tsv", "--out-summary", tmp_path / "ts.tsv")
    t = {r["seq_id"]: r for r in table(tmp_path / "t.tsv")}
    assert (t["chr1"]["telomere_start"], t["chr1"]["telomere_end"]) == ("yes", "yes")
    assert (t["chr2"]["telomere_start"], t["chr2"]["telomere_end"]) == ("yes", "no")
    assert (t["chr3"]["telomere_start"], t["chr3"]["telomere_end"]) == ("no", "no")
    s = dict(l.split("\t") for l in (tmp_path / "ts.tsv").read_text().split("\n") if l)
    assert s["telomere_both_ends"] == "1" and s["telomere_one_end"] == "1"


def test_summary_marks_missing_as_not_run_and_parses_present(tmp_path):
    (tmp_path / "s.tsv").write_text("metric\tscaffold\tcontig\nN50\t100\t50\ngaps\t1\t-\n")
    (tmp_path / "b.txt").write_text("\tC:97.0%[S:95.0%,D:2.0%],F:1.0%,M:2.0%,n:2326\n")
    (tmp_path / "x.LAI").write_text("Chr\tFrom\tTo\tIntact\tTotal\traw_LAI\tLAI\nwhole_genome\t1\t1000\t0.05\t0.4\t12.5\t14.07\n")
    (tmp_path / "g.qv").write_text("genome\t1000\t5000000\t42.3\t5.9e-05\n")
    (tmp_path / "g.cmp").write_text("genome\tall\t900\t1000\t96.5\n")
    (tmp_path / "fcs.txt").write_text("##[...]\n#seq_id\tstart_pos\tend_pos\tseq_len\taction\n"
                                      "s9\t1\t500\t500\tEXCLUDE\ts1\t1\t50\t500\tTRIM\n".replace("\ts1", "\ns1"))
    run("assembly_qc_summary.py", "--stats", tmp_path / "s.tsv", "--busco", f"eud={tmp_path / 'b.txt'}",
        "--lai", tmp_path / "x.LAI", "--merqury-qv", tmp_path / "g.qv", "--merqury-completeness", tmp_path / "g.cmp",
        "--fcs", tmp_path / "fcs.txt", "--out", tmp_path / "o.tsv")
    m = {r["metric"]: r["value"] for r in table(tmp_path / "o.tsv")}
    assert m["scaffold_N50"] == "100" and m["contig_N50"] == "50" and m["gaps"] == "1"
    assert m["BUSCO_genome_eud"].startswith("C:97.0%")
    assert m["LAI"] == "raw_LAI=12.5, LAI=14.07"
    assert m["Merqury_QV"].startswith("QV=42.3") and m["Merqury_kmer_completeness"] == "96.5%"
    assert m["tidk"] == "not run"
    assert m["FCS-GX_flagged_regions"] == "EXCLUDE=1, TRIM=1"
    run("assembly_qc_summary.py", "--stats", tmp_path / "s.tsv", "--out", tmp_path / "o2.tsv")
    m = {r["metric"]: r["value"] for r in table(tmp_path / "o2.tsv")}
    assert m["Merqury_QV"] == m["FCS-GX"] == "not run"
    assert "LAI" not in m                                    # LAI is only listed when it was run


def test_toggles(tmp_path):
    p = make(tmp_path)
    # defaults: telomeres on with TTTAGGG, LAI on (after EDTA), Merqury/FCS-GX off without data
    assert p.run_telomere and p.telomere_motif == "TTTAGGG" and p.run_lai
    assert not p.run_merqury and not p.run_fcs_gx
    assert not make(tmp_path, assembly_qc={"lai": False}).run_lai
    p = make(tmp_path, assembly_qc={"lai": True, "merqury": {"reads": ["wgs.fq.gz"]},
                                    "telomere": {"motif": "ttaggg"}, "fcs_gx": {"db": "gx/all", "tax_id": 3702}})
    assert p.run_lai and p.run_merqury and p.run_fcs_gx
    assert p.telomere_motif == "TTAGGG"                      # user-specified sequence
    assert not make(tmp_path, assembly_qc={"telomere": {"enabled": False}}).run_telomere


def test_toggle_validation(tmp_path):
    with pytest.raises(ConfigError, match="fcs_gx"):
        make(tmp_path, assembly_qc={"fcs_gx": {"db": "gx/all"}})
    with pytest.raises(ConfigError, match="motif"):
        make(tmp_path, assembly_qc={"telomere": {"motif": "TTAGGGX"}})


def test_summary_lai_not_applicable(tmp_path):
    (tmp_path / "s.tsv").write_text("metric\tscaffold\tcontig\nN50\t100\t50\n")
    (tmp_path / "x.LAI").write_text("# LAI not applicable: Total LTR sequence content (1.92%) is too low for "
                                    "accurate LAI calculation (min 5% required)\n")
    run("assembly_qc_summary.py", "--stats", tmp_path / "s.tsv", "--lai", tmp_path / "x.LAI", "--out", tmp_path / "o.tsv")
    m = {r["metric"]: r["value"] for r in table(tmp_path / "o.tsv")}
    assert m["LAI"] == "not applicable (LTR content 1.92% < 5%)"
