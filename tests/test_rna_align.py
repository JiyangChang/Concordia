"""rna_align_stats.py on realistic HISAT2 / samtools output."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "rna_align_stats.py")

HISAT2_PE = """10000 reads; of these:
  10000 (100.00%) were paired; of these:
    1000 (10.00%) aligned concordantly 0 times
    8500 (85.00%) aligned concordantly exactly 1 time
    500 (5.00%) aligned concordantly >1 times
    ----
    1000 pairs aligned concordantly 0 times; of these:
      100 (10.00%) aligned discordantly 1 time
95.00% overall alignment rate
"""
HISAT2_SE = """5000 reads; of these:
  5000 (100.00%) were unpaired; of these:
    500 (10.00%) aligned 0 times
    4000 (80.00%) aligned exactly 1 time
    500 (10.00%) aligned >1 times
90.00% overall alignment rate
"""
FLAGSTAT = """20000 + 0 in total (QC-passed reads + QC-failed reads)
19000 + 0 primary
1000 + 0 secondary
0 + 0 supplementary
18100 + 0 mapped (90.50% : N/A)
17100 + 0 primary mapped (90.00% : N/A)
"""


def sn(total, mapped, err, length):
    return (f"raw total sequences:\t{total}\nreads mapped:\t{mapped}\n"
            f"error rate:\t{err}\t# mismatches / bases mapped (cigar)\naverage length:\t{length}\n")


def test_rna_align_stats(tmp_path):
    q = tmp_path / "qc"
    for k in ("short", "long"):
        (q / k).mkdir(parents=True)
    (q / "short" / "leaf.samtools_stats.txt").write_text(sn(20000, 19000, "2.1e-03", 150))
    (q / "short" / "leaf.spliced.txt").write_text("4750\n")
    (q / "leaf.hisat2.txt").write_text(HISAT2_PE)
    (q / "short" / "root.samtools_stats.txt").write_text(sn(5000, 4500, "2.5e-03", 100))
    (q / "short" / "root.spliced.txt").write_text("900\n")
    (q / "root.hisat2.txt").write_text(HISAT2_SE)
    (q / "long" / "iso.samtools_stats.txt").write_text(sn(1000, 980, "1.2e-02", 2500))
    (q / "long" / "iso.spliced.txt").write_text("882\n")
    (q / "sr.flagstat.txt").write_text(FLAGSTAT)
    out_s, out_m = tmp_path / "s.tsv", tmp_path / "m.tsv"
    r = subprocess.run([sys.executable, SCRIPT, "--qc-dir", q,
                        "--samples", "short:leaf:illumina", "short:root:illumina", "long:iso:pacbio_isoseq",
                        "--merged", "sr=BRAKER4 (bam_files);GeMoMa ERE",
                        "--out-samples", out_s, "--out-merged", out_m], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = {x["sample"]: x for x in csv.DictReader(open(out_s), delimiter="\t")}
    assert rows["leaf"]["mapping_%"] == "95.00" and rows["leaf"]["unique_%"] == "85.00" and rows["leaf"]["multi_%"] == "5.00"
    assert rows["leaf"]["spliced_%"] == "25.00"
    assert rows["root"]["unique_%"] == "80.00" and rows["root"]["multi_%"] == "10.00"     # single-end summary
    assert rows["iso"]["type"] == "long" and rows["iso"]["unique_%"] == "-" and rows["iso"]["spliced_%"] == "90.00"
    assert rows["iso"]["error_rate"] == "1.2e-02" and rows["iso"]["mean_length"] == "2500"
    m = list(csv.DictReader(open(out_m), delimiter="\t"))[0]
    assert m["bam"] == "sr.bam" and m["primary_mapped"] == "17100" and m["primary_mapped_%"] == "90.00"
    assert m["used_by"] == "BRAKER4 (bam_files); GeMoMa ERE"
