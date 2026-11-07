"""function_stats.py on realistic InterProScan / eggNOG-mapper output."""
import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "workflow", "scripts", "function_stats.py")

GFF = "".join(f"c1\tx\tgene\t{i * 1000}\t{i * 1000 + 500}\t.\t+\t.\tID=G{i}\n"
              f"c1\tx\tmRNA\t{i * 1000}\t{i * 1000 + 500}\t.\t+\t.\tID=G{i}.t1;Parent=G{i}\n" for i in range(1, 6))


def ipr(pid, db, acc, desc, ipr_acc="-", ipr_desc="-", go="-"):
    return "\t".join([pid, "md5", "300", db, acc, desc, "1", "100", "1e-10", "T", "01-01-2026", ipr_acc, ipr_desc, go]) + "\n"


IPR = (ipr("G1.t1", "Pfam", "PF00069", "Protein kinase domain", "IPR000719", "Protein kinase domain",
           "GO:0004672(InterPro)|GO:0006468(InterPro)") +
       ipr("G1.t1", "PANTHER", "PTHR1", "x") +
       ipr("G2.t1", "Pfam", "PF00069", "Protein kinase domain", "IPR000719", "Protein kinase domain") +
       ipr("G3.t1", "Gene3D", "G3D.1", "x"))
EGG_HEADER = "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\n"
EGG = ("## emapper-2.1\n" + EGG_HEADER +
       "G1.t1\ts\t1e-50\t200\tOG\t33090|Viridiplantae\tT\tkinase\tKIN1\tGO:0004672\t2.7.11.1\tko:K04730\tko04010,map04010\n"
       "G2.t1\ts\t1e-50\t200\tOG\t33090|Viridiplantae\tTK\tkinase\t-\t-\t-\t-\t-\n"
       "G4.t1\ts\t1e-50\t200\tOG\t2759|Eukaryota\tS\tunknown\t-\t-\t-\tko:K00001\tmap00010\n")
OBO = """[Term]
id: GO:0004672
name: protein kinase activity
namespace: molecular_function

[Term]
id: GO:0006468
name: protein phosphorylation
namespace: biological_process
"""


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def run(tmp_path, obo=True, kegg=True):
    (tmp_path / "f.gff3").write_text(GFF)
    (tmp_path / "c1.tsv").write_text(IPR)
    (tmp_path / "e.tsv").write_text(EGG)
    (tmp_path / "go.obo").write_text(OBO)
    out = tmp_path / "stats"
    args = [sys.executable, SCRIPT, "--gff", tmp_path / "f.gff3", "--interproscan", tmp_path / "c1.tsv",
            tmp_path / "missing_chunk.tsv", "--eggnog", tmp_path / "e.tsv", "--outdir", out]
    if obo:
        args += ["--go-obo", tmp_path / "go.obo"]
    if kegg:
        (tmp_path / "kegg.tsv").write_text("map04010\tMAPK signaling pathway\nmap00010\tGlycolysis / Gluconeogenesis\n")
        args += ["--kegg-names", tmp_path / "kegg.tsv"]
    r = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return out


def test_function_stats(tmp_path):
    out = run(tmp_path)
    ov = {x["category"]: x for x in rows(out / "overview.tsv")}
    assert ov["total genes"]["genes"] == "5"
    assert ov["InterProScan (any member database)"]["genes"] == "3"           # G1 G2 G3
    assert ov["eggNOG-mapper"]["genes"] == "3"                                # G1 G2 G4
    assert ov["both"]["genes"] == "2" and ov["InterProScan only"]["genes"] == "1"
    assert ov["eggNOG-mapper only"]["genes"] == "1" and ov["neither (unannotated)"]["genes"] == "1"   # G5
    assert ov["either (annotated)"]["percent"] == "80.00" and ov["with GO term(s)"]["genes"] == "1"
    db = {x["database"]: x["genes"] for x in rows(out / "interpro_databases.tsv")}
    assert db == {"Pfam": "2", "PANTHER": "1", "Gene3D": "1"}
    pf = rows(out / "top_pfam.tsv")[0]
    assert (pf["pfam"], pf["description"], pf["genes"]) == ("PF00069", "Protein kinase domain", "2")
    go = {x["go"]: x for x in rows(out / "go_top.tsv")}
    assert go["GO:0004672"]["name"] == "protein kinase activity" and go["GO:0004672"]["namespace"] == "MF"
    ns = {x["namespace"]: x["genes"] for x in rows(out / "go_namespace.tsv")}
    assert ns == {"BP": "1", "MF": "1", "CC": "0"}
    kg = {x["category"]: x["genes"] for x in rows(out / "kegg.tsv")}
    assert kg == {"with KEGG KO": "2", "with KEGG pathway": "2", "with EC number": "1"}
    assert {x["pathway"]: (x["name"], x["genes"]) for x in rows(out / "kegg_top_pathways.tsv")} == \
        {"map04010": ("MAPK signaling pathway", "1"), "map00010": ("Glycolysis / Gluconeogenesis", "1")}
    cog = {x["code"]: x for x in rows(out / "cog_categories.tsv")}
    assert cog["T"]["genes"] == "2" and cog["T"]["description"] == "Signal transduction mechanisms"
    assert cog["K"]["genes"] == "1" and cog["S"]["genes"] == "1"
    assert {x["max_annot_lvl"]: x["genes"] for x in rows(out / "eggnog_levels.tsv")} == \
        {"33090|Viridiplantae": "2", "2759|Eukaryota": "1"}


def test_function_stats_without_obo(tmp_path):
    out = run(tmp_path, obo=False, kegg=False)
    assert {x["name"] for x in rows(out / "kegg_top_pathways.tsv")} == {"-"}
    assert rows(out / "go_namespace.tsv") == []
    assert rows(out / "go_top.tsv")[0]["name"] == "-"
