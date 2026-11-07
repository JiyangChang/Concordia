"""reference_compare.py: CDS-level agreement with a reference annotation."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "reference_compare.py")


def g(seq, gid, cds, strand="+"):
    rows = [f"{seq}\tx\tgene\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\tID={gid}",
            f"{seq}\tx\tmRNA\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\tID={gid}.t;Parent={gid}"]
    rows += [f"{seq}\tx\tCDS\t{s}\t{e}\t.\t{strand}\t0\tParent={gid}.t" for s, e in cds]
    return "\n".join(rows) + "\n"


def test_reference_compare(tmp_path):
    # reference uses original names (Chr1, Chr2); the analysed genome only contains Chr1 -> seq0001
    ref = (g("Chr1", "r1", [(100, 200), (300, 400)]) + g("Chr1", "r2", [(1000, 1300)]) +
           g("Chr1", "r3", [(5000, 5100), (5200, 5300)]) + g("Chr2", "rX", [(100, 400)]))
    pred = (g("seq0001", "p1", [(100, 200), (300, 400)]) +          # exact
            g("seq0001", "p3", [(5000, 5100), (5250, 5300)]) +      # overlaps r3, different intron
            g("seq0001", "p9", [(9000, 9300)]))                     # novel
    (tmp_path / "ref.gff3").write_text(ref)
    (tmp_path / "pred.gff3").write_text(pred)
    (tmp_path / "map.tsv").write_text("new_id\toriginal_id\nseq0001\tChr1\n")
    out = tmp_path / "o.tsv"
    r = subprocess.run([sys.executable, SCRIPT, "--reference", tmp_path / "ref.gff3", "--seqid-map", tmp_path / "map.tsv",
                        "--set", f"final={tmp_path / 'pred.gff3'}", "--out", out], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    row = list(csv.DictReader((l for l in open(out) if not l.startswith("#")), delimiter="\t"))[0]
    assert row["ref_genes"] == "3"                                   # Chr2 not analysed
    assert row["gene_recall_%"] == "66.7" and row["gene_precision_%"] == "66.7"
    assert row["intron_recall_%"] == "50.0" and row["intron_precision_%"] == "50.0"
    assert row["tx_recall_%"] == "33.3" and row["tx_precision_%"] == "33.3"
