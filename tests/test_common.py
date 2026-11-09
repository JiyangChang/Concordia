import os
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


@pytest.mark.parametrize("short,pacbio,ont,mode,tags", [
    (True, True, True, "dual", ["bam", "bam_isoseq", "bam_isoseq"]),
    (True, True, False, "dual", ["bam", "bam_isoseq"]),
    (True, False, True, "ETP", ["bam", "bam_isoseq"]),      # ONT never goes to BRAKER4
    (True, False, False, "ETP", ["bam"]),
    (False, True, False, "IsoSeq", ["bam_isoseq"]),
    (False, False, True, "EP", ["bam_isoseq"]),             # only-ONT: provisional bam_isoseq
    (False, True, True, "IsoSeq", ["bam_isoseq", "bam_isoseq"]),
])
def test_modes(tmp_path, short, pacbio, ont, mode, tags):
    p = make(tmp_path, short, pacbio, ont)
    assert p.braker4_mode() == mode
    lines, _ = p.eviann_reads()
    assert [t for _, t in lines] == tags                    # one EviAnn line per sample, per-sample BAMs
    assert all("/02_RNA/align_" in b and "/merged/" not in b for bams, _ in lines for b in bams)
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
    assert u["sr"] == ["BRAKER4 (bam_files)", "GeMoMa ERE"]          # EviAnn uses the per-sample BAMs
    assert u["pacbio"] == ["BRAKER4 (isoseq_bam)"]
    assert u["lr"] == ["StringTie -> EVM TRANSCRIPT (LR_stringtie)"]
    assert u["ont"] == ["merged into lr.bam (not used by BRAKER4)"]
    u = make(tmp_path, short=True).bam_usage()
    assert u == {"sr": ["BRAKER4 (bam_files)", "StringTie -> EVM TRANSCRIPT (SR_stringtie)", "GeMoMa ERE"]}
    u = make(tmp_path, short=False, ont=True).bam_usage()
    assert u["lr"] == ["StringTie -> EVM TRANSCRIPT (LR_stringtie)"]


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


def test_path_settings_tilde_and_relative(tmp_path, monkeypatch):
    import os
    from common import _abs
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    base = str(tmp_path / "analysis")
    assert _abs("~/db/eggnog", base) == str(tmp_path / "home" / "db" / "eggnog")
    assert _abs("db/dfam", base) == os.path.join(base, "db", "dfam")
    assert _abs("/abs/x", base) == "/abs/x"
    for unset in (None, "", "null"):
        assert _abs(unset, base) is None


def test_resource_units():
    import pytest
    from common import ConfigError, normalize_resources, parse_mem_mb, parse_minutes
    # binary units as sbatch --mem: 16G = 16384 MB
    assert parse_mem_mb("16G") == 16384 and parse_mem_mb("16GB") == 16384 and parse_mem_mb("16GiB") == 16384
    assert parse_mem_mb("512M") == 512 and parse_mem_mb("1.5T") == 1572864 and parse_mem_mb("0.5g") == 512
    assert parse_mem_mb(32000) == 32000 and parse_mem_mb(16, "gb") == 16384
    assert parse_minutes(720) == 720 and parse_minutes("12h") == 720 and parse_minutes("3d") == 4320
    assert parse_minutes("1d12h") == 2160 and parse_minutes("90m") == 90 and parse_minutes("2h30min") == 150
    assert normalize_resources({"threads": 8, "mem": "16G", "runtime": "12h"}) == {"threads": 8, "mem_mb": 16384, "runtime": 720}
    assert normalize_resources({"mem": 64}) == {"mem_mb": 65536}                  # bare number in mem = G
    assert normalize_resources({"mem_mb": 8000, "runtime": 600}) == {"mem_mb": 8000, "runtime": 600}   # old style
    for bad in ({"mem": "16X"}, {"runtime": "12 hours"}, {"mem": "16G", "mem_mb": 16000}):
        with pytest.raises(ConfigError):
            normalize_resources(bad)


def test_resources_step_over_default(tmp_path):
    from common import Project
    cfg = {"resources": {"default": {"threads": 4, "mem": "8G", "runtime": "10h"}, "edta": {"threads": 48, "mem": "128G"}}}
    p = Project.__new__(Project)
    p.cfg = cfg
    assert p.resources("edta") == {"threads": 48, "mem_mb": 131072, "runtime": 600}
    assert p.resources("other") == {"threads": 4, "mem_mb": 8192, "runtime": 600}


def test_eviann_tissue_groups(tmp_path):
    from common import ConfigError
    s = ("sample\tr1\tr2\tstrandness\ttissue\n"
         "leaf1\ta\tb\tRF\tleaf\nleaf2\ta\tb\tRF\tleaf\nroot1\ta\tb\tRF\troot\nfree\ta\tb\tRF\t\n")
    l = "sample\treads\tplatform\tsample_group\nisoL\tx.fa\tpacbio_isoseq\tleaf\nisoX\ty.fa\tpacbio_isoseq\t\n"
    p = Project({"genome": "g.fa", "outdir": "out", "homologs": write(tmp_path, "h.tsv", "name\tgenome\tgff3\n"
                                                                    "ath\tg\tf\n"),
                 "rna": {"short": write(tmp_path, "s.tsv", s), "long": write(tmp_path, "l.tsv", l)},
                 "qc": {"busco_lineages": ["x"], "omark_db": "LUCA.h5"}}, str(tmp_path))
    lines, summary = p.eviann_reads()
    assert [(len(b), t) for b, t in lines] == [(1, "bam"), (1, "bam_isoseq"), (2, "bam_mix"), (1, "bam")]
    leaf = next(b for b, t in lines if t == "bam_mix")
    assert [os.path.relpath(b, p.out) for b in leaf] == ["02_RNA/EviAnn/groups/leaf.sr.bam",   # 2 samples: merged
                                                         "02_RNA/align_long/isoL.bam"]          # 1 sample: as is
    root = lines[-1][0][0]
    assert os.path.relpath(root, p.out) == "02_RNA/align_short/root1.bam"
    assert summary == "2 x bam + 1 x bam_isoseq + 1 x bam_mix"
    assert p.eviann_groups()["leaf"]["sr"] == [p.o("02_RNA/align_short/leaf1.bam"), p.o("02_RNA/align_short/leaf2.bam")]
    with pytest.raises(ConfigError, match="tissue"):
        Project({"genome": "g.fa", "rna": {"short": write(tmp_path, "b.tsv", "sample\tr1\tr2\tstrandness\ttissue\n"
                                                                             "a\tx\ty\tRF\tleaf tissue\n")}},
                str(tmp_path))
