"""BRAKER4 gene_support summary and GeMoMa evidence statistics."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "workflow", "scripts")


def run(name, *args):
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, name), *map(str, args)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


def rows(path):
    return list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))


def test_braker4_support(tmp_path):
    (tmp_path / "s.tsv").write_text("gene_id\tRNA_support\tprotein_support\tsource\n"
                                    "g1\t1.0\t0.5\tyes\ng2\t0\t1\tno\ng3\t0.25\t0\tyes\ng4\t0\t0\tno\n")
    out = tmp_path / "o.tsv"
    run("braker4_support.py", "--support", tmp_path / "s.tsv", "--mode", "dual", "--out", out)
    r = {x["column"]: x for x in rows(out)}
    assert r["RNA_support"]["supported"] == "2" and r["RNA_support"]["supported_%"] == "50.0"
    assert r["RNA_support"]["mean"] == "0.312"
    assert r["protein_support"]["supported_%"] == "50.0"
    assert r["source"]["supported"] == "2" and r["source"]["mean"] == "-"          # text column
    assert "# BRAKER4 mode dual; 4 models" in out.read_text()


def test_braker4_support_missing(tmp_path):
    (tmp_path / "e.tsv").write_text("")
    out = tmp_path / "o.tsv"
    run("braker4_support.py", "--support", tmp_path / "e.tsv", "--mode", "ETP", "--out", out)
    assert "missing or empty" in out.read_text() and rows(out) == []


def test_gemoma_evidence(tmp_path):
    mrna = "c1\tGeMoMa\tmRNA\t1\t100\t.\t+\t.\tID={id};evidence={ev};tie={tie};tpc={tpc};score={sc}{alt}\n"
    final = tmp_path / "merged.gff"
    final.write_text(mrna.format(id="spA_t1", ev=2, tie=1, tpc=1, sc=300, alt=";alternative=\"spB_t9\"") +
                     mrna.format(id="spA_t2", ev=1, tie=0.5, tpc=0.8, sc=200, alt="") +
                     mrna.format(id="spB_t3", ev=1, tie=0, tpc=0, sc=100, alt=""))
    pred = "c1\tGeMoMa\tmRNA\t1\t100\t.\t+\t.\tID=x{}\n"
    files = {}
    for sp, n_unf, n_fin in (("spA", 3, 2), ("spB", 1, 1)):
        u, f = tmp_path / f"{sp}.unf.gff", tmp_path / f"{sp}.final.gff"
        u.write_text("".join(pred.format(i) for i in range(n_unf)))
        f.write_text("".join(pred.format(i) for i in range(n_fin)))
        files[sp] = f"{sp}={u},{f}"
    refs, attrs = tmp_path / "r.tsv", tmp_path / "a.tsv"
    run("gemoma_evidence.py", "--final", final, "--reference", files["spA"], "--reference", files["spB"],
        "--out-references", refs, "--out-attributes", attrs)
    r = {x["species"]: x for x in rows(refs)}
    assert r["spA"]["unfiltered_predictions"] == "3" and r["spB"]["unfiltered_predictions"] == "1"
    assert r["spA"]["species_final_transcripts"] == "2" and r["spB"]["species_final_transcripts"] == "1"
    assert r["spA"]["merged_representative"] == "2" and r["spB"]["merged_representative"] == "1"
    assert r["spA"]["merged_supported"] == "2" and r["spB"]["merged_supported"] == "2"   # spB via alternative
    m = {x["metric"]: x["value"] for x in rows(attrs)}
    assert m["merged_transcripts"] == "3"
    assert m["transcripts_supported_by_1_species"] == "2" and m["transcripts_supported_by_2_species"] == "1"
    assert m["tie (intron RNA support) mean"] == "0.500"
    assert m["tie (intron RNA support) = 1 (fully supported) %"] == "33.3"
    assert m["tpc (RNA coverage) = 0 (no support) %"] == "33.3"


def test_gemoma_evidence_without_rna_attributes(tmp_path):
    final = tmp_path / "f.gff"
    final.write_text("c1\tGeMoMa\tmRNA\t1\t100\t.\t+\t.\tID=t1;score=50\n")
    refs, attrs = tmp_path / "r.tsv", tmp_path / "a.tsv"
    run("gemoma_evidence.py", "--final", final, "--reference", f"spA={tmp_path / 'missing.gff'},",
        "--out-references", refs, "--out-attributes", attrs)
    m = {x["metric"]: x["value"] for x in rows(attrs)}
    assert m["tie (intron RNA support)"] == "not available"
    assert rows(refs)[0]["unfiltered_predictions"] == "-"


def gff_tx(gid, cds, strand="+", seq="c1"):
    lo, hi = cds[0][0], cds[-1][1]
    out = [f"{seq}\tx\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID={gid}",
           f"{seq}\tx\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID={gid}.t;Parent={gid}"]
    out += [f"{seq}\tx\tCDS\t{s}\t{e}\t.\t{strand}\t0\tParent={gid}.t" for s, e in cds]
    return "\n".join(out) + "\n"


def test_source_support(tmp_path):
    q = tmp_path / "q.gff3"
    q.write_text(gff_tx("e1", [(100, 200), (300, 400)]) +      # exact in A and B
                 gff_tx("e2", [(1000, 1300)]) +                 # exact in A, overlaps C
                 gff_tx("e3", [(5000, 5300)], strand="-") +     # overlaps A only on the other strand -> none
                 gff_tx("e4", [(9000, 9300)]))                  # nothing
    (tmp_path / "a.gff3").write_text(gff_tx("a1", [(100, 200), (300, 400)]) + gff_tx("a2", [(1000, 1300)]) +
                                     gff_tx("a3", [(5000, 5300)], strand="+"))
    (tmp_path / "b.gff3").write_text(gff_tx("b1", [(100, 200), (300, 400)]))
    (tmp_path / "c.gff3").write_text(gff_tx("c1", [(1200, 1500)]))
    sup, comb = tmp_path / "s.tsv", tmp_path / "c.tsv"
    run("source_support.py", "--query", q, "--source", f"A={tmp_path / 'a.gff3'}", "--source", f"B={tmp_path / 'b.gff3'}",
        "--source", f"C={tmp_path / 'c.gff3'}", "--out-support", sup, "--out-combos", comb)
    s = {x["source"]: x for x in rows(sup)}
    assert (s["A"]["exact"], s["A"]["overlap"]) == ("2", "2")
    assert (s["B"]["exact"], s["C"]["exact"], s["C"]["overlap"]) == ("1", "0", "1")
    assert s["none (no source overlaps)"]["overlap"] == "2" and s["A"]["exact_%"] == "50.0"
    c = {x["combination"]: x["transcripts"] for x in rows(comb)}
    assert c == {"A+B": "1", "A": "1", "no exact match": "2"}
