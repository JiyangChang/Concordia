"""Decision logic in workflow/lib/common.py."""
import pytest

from common import ConfigError, Project


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def make(tmp_path, short=True, pacbio=False, ont=False, **extra):
    s = "sample\tr1\tr2\tstrandness\n" + ("leaf\ta_R1.fq\ta_R2.fq\tRF\n" if short else "")
    l = "sample\treads\tplatform\n"
    if pacbio:
        l += "iso\tiso.fq\tpacbio_isoseq\n"
    if ont:
        l += "ont\tont.fq\tont_cdna\n"
    cfg = {
        "outdir": "out", "genome": "g.fa", "species_prefix": "AB",
        "rna": {"short": write(tmp_path, "s.tsv", s), "long": write(tmp_path, "l.tsv", l)},
        "homologs": write(tmp_path, "h.tsv", "name\tgenome\tgff3\tprotein\nspA\ta.fa\ta.gff3\t\n"),
        "qc": {"busco_lineages": ["eudicots_odb12"], "omark_db": "LUCA.h5"},
        "evm": {"weights": {"ABINITIO_PREDICTION": {"BRAKER4": 2},
                            "TRANSCRIPT": {"LR_stringtie": 10, "SR_stringtie": 6}}},
    }
    cfg.update(extra)
    return Project(cfg, str(tmp_path))


@pytest.mark.parametrize("short,pacbio,ont,mode,tag", [
    (True, True, True, "dual", "bam_mix"),
    (True, True, False, "dual", "bam_mix"),
    (True, False, True, "ETP", "bam_mix"),      # ONT never goes to BRAKER4
    (True, False, False, "ETP", "bam"),
    (False, True, False, "IsoSeq", "bam_isoseq"),
    (False, False, True, "EP", "bam_isoseq"),   # only-ONT: provisional bam_isoseq
    (False, True, True, "IsoSeq", "bam_isoseq"),
])
def test_modes(tmp_path, short, pacbio, ont, mode, tag):
    p = make(tmp_path, short, pacbio, ont)
    assert p.braker4_mode() == mode
    assert p.eviann_reads()[1] == tag
    bam, iso = p.braker4_bams()
    assert bool(iso) == pacbio and (not iso or iso.endswith("pacbio.bam"))
    assert bool(bam) == short


def test_transcript_evidence_and_weights(tmp_path):
    p = make(tmp_path, short=True)
    assert p.transcript_evidence()[0] == "SR_stringtie"
    assert "TRANSCRIPT\tSR_stringtie\t6" in p.evm_weights()
    assert not any("LR_stringtie" in w for w in p.evm_weights())
    p = make(tmp_path, short=True, ont=True)
    assert p.transcript_evidence() == ("LR_stringtie", "lr", True)


def test_busco_lineage_required(tmp_path):
    with pytest.raises(ConfigError, match="busco_lineages"):
        make(tmp_path, qc={"busco_lineages": [], "omark_db": "x"})


def test_no_rna_rejected(tmp_path):
    with pytest.raises(ConfigError, match="rna.allow_none"):     # error names the opt-in switch
        make(tmp_path, short=False)


def test_no_rna_allowed(tmp_path):
    p = make(tmp_path, short=False, rna={"short": None, "long": None, "allow_none": True},
             evm={"weights": {"ABINITIO_PREDICTION": {"BRAKER4": 2}, "OTHER_PREDICTION": {"GeMoMa": 5, "EviAnn": 8},
                              "PROTEIN": {"miniprot": 3}, "TRANSCRIPT": {"LR_stringtie": 10, "SR_stringtie": 6}},
                  "tune": {"grid": {"BRAKER4": [1, 2], "EviAnn": [5, 10], "transcripts": [5, 10]}}})
    assert not p.has_rna and p.braker4_mode() == "EP" and p.braker4_bams() == ("", "")
    assert p.transcript_evidence() == (None, None, False) and p.eviann_reads() == ([], "none")
    assert p.gemoma_ere_bam() is None and p.bam_usage() == {}
    assert p.evm_weights() == ["ABINITIO_PREDICTION\tBRAKER4\t2", "OTHER_PREDICTION\tGeMoMa\t5", "PROTEIN\tminiprot\t3"]
    assert [c[1] for c in p.evm_grid()] == [{"BRAKER4": 2, "GeMoMa": 5, "miniprot": 3},
                                            {"BRAKER4": 1, "GeMoMa": 5, "miniprot": 3}]


def test_bad_platform(tmp_path):
    l = write(tmp_path, "bad.tsv", "sample\treads\tplatform\nx\tx.fq\tnanopore\n")
    with pytest.raises(ConfigError, match="platform"):
        make(tmp_path, rna={"short": write(tmp_path, "s2.tsv", "sample\tr1\tr2\tstrandness\n"), "long": l})


def test_edta_u(tmp_path):
    assert make(tmp_path).edta_u_arg() == ""
    assert make(tmp_path, repeat={"edta_u": 9.12e-9}).edta_u_arg() == "-u 9.12e-09"


def test_evm_overlap(tmp_path):
    p = make(tmp_path, max_intron=50000)
    assert p.evm_overlap() == (500000, 100000)
    p = make(tmp_path, max_intron=500000)
    assert p.evm_overlap() == (500000, 250000)   # capped below segment size


def test_gemoma_ere(tmp_path):
    assert make(tmp_path).gemoma_ere_bam().endswith("sr.bam")
    assert make(tmp_path, short=False, pacbio=True).gemoma_ere_bam() is None
    assert make(tmp_path, short=False, pacbio=True, gemoma={"ere_long_reads": True}).gemoma_ere_bam().endswith("lr.bam")


def test_bam_usage(tmp_path):
    u = make(tmp_path, short=True, pacbio=True, ont=True).bam_usage()
    assert u["sr"] == ["BRAKER4 (bam_files)", "EviAnn (bam_mix)", "GeMoMa ERE"]
    assert u["pacbio"] == ["BRAKER4 (isoseq_bam)"]
    assert u["lr"] == ["EviAnn (bam_mix)", "StringTie -> EVM TRANSCRIPT (LR_stringtie)"]
    assert u["ont"] == ["merged into lr.bam (not used by BRAKER4)"]
    u = make(tmp_path, short=True).bam_usage()
    assert u == {"sr": ["BRAKER4 (bam_files)", "EviAnn (bam)", "StringTie -> EVM TRANSCRIPT (SR_stringtie)", "GeMoMa ERE"]}
    u = make(tmp_path, short=False, ont=True).bam_usage()
    assert u["lr"] == ["EviAnn (bam_isoseq)", "StringTie -> EVM TRANSCRIPT (LR_stringtie)"]


def test_evm_grid(tmp_path):
    p = make(tmp_path, short=True, pacbio=True, evm={"weights": {
        "ABINITIO_PREDICTION": {"BRAKER4": 2}, "OTHER_PREDICTION": {"GeMoMa": 5, "EviAnn": 8},
        "PROTEIN": {"miniprot": 3}, "TRANSCRIPT": {"LR_stringtie": 10, "SR_stringtie": 6}}})
    g = p.evm_grid()
    assert g[0] == ("c00", {"BRAKER4": 2, "GeMoMa": 5, "EviAnn": 8, "miniprot": 3, "LR_stringtie": 10})
    assert len(g) == 1 + 2 * 2 * 2 * 1 * 2                      # default grid, baseline not duplicated
    assert all(c[1]["miniprot"] == 2 for c in g[1:]) and {c[1]["LR_stringtie"] for c in g[1:]} == {5, 10}
    assert p.evm_weight_lines(g[1][1])[0] == "ABINITIO_PREDICTION\tBRAKER4\t1"
    # a grid containing the baseline values does not repeat it
    p = make(tmp_path, short=True, evm={"weights": {"ABINITIO_PREDICTION": {"BRAKER4": 2},
                                                    "TRANSCRIPT": {"SR_stringtie": 6}},
                                        "tune": {"grid": {"BRAKER4": [1, 2], "transcripts": [6]}}})
    assert [c[1] for c in p.evm_grid()] == [{"BRAKER4": 2, "SR_stringtie": 6}, {"BRAKER4": 1, "SR_stringtie": 6}]
    assert make(tmp_path, evm={"tune": {"enabled": False}}).evm_tune["enabled"] is False


def test_long_read_counts_column(tmp_path):
    l = write(tmp_path, "long_counts.tsv", "sample\treads\tplatform\tcounts\niso\tiso.fq\tpacbio_isoseq\tiso.cluster_report.csv\n"
                                  "ont\tont.fq\tont_cdna\t\n")
    p = make(tmp_path, rna={"short": write(tmp_path, "s.tsv", "sample\tr1\tr2\tstrandness\n"), "long": l})
    assert p.long_sample("iso")["counts"] == str(tmp_path / "iso.cluster_report.csv")
    assert p.long_sample("ont")["counts"] is None
