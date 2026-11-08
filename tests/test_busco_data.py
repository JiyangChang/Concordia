"""BUSCO dataset directory choice (lib/busco_data.py) and the dataset fetcher (scripts/fetch_busco.py)."""
import os
import stat

import pytest

import busco_data
import fetch_busco

VERSIONS = ("brassicales_odb12\t2026-05-22\tabc\tEukaryota\tlineages\n"
            "eudicotyledons_odb12\t2026-05-22\tdef\tEukaryota\tlineages\n"
            "eukaryota_odb12\t2026-05-22\t123\tEukaryota\tlineages\n"
            "list_of_reference_markers.eukaryota_odb12.txt\t2025-01-15\t456\tEukaryota\tplacement_files\n"
            "tree.eukaryota_odb12.2.nwk\t2026-05-13\t789\tEukaryota\tplacement_files\n"
            "tree.bacteria_odb12.nwk\t2025-01-15\t000\tProkaryota\tplacement_files\n")


def dataset(d, lin):
    os.makedirs(d / "lineages" / lin, exist_ok=True)
    (d / "lineages" / lin / "dataset.cfg").write_text("name=" + lin + "\n")


@pytest.fixture
def readonly(tmp_path):
    d = tmp_path / "shared"
    dataset(d, "brassicales_odb12")
    os.chmod(d, stat.S_IRUSR | stat.S_IXUSR)
    yield d
    os.chmod(d, stat.S_IRWXU)


def test_directory_choice(tmp_path, readonly):
    local = str(tmp_path / "out" / "busco_downloads")
    assert busco_data.writable(local)                       # does not exist yet: parent decides
    assert busco_data.lineage_dir(None, local, "brassicales_odb12") == local
    if os.access(readonly, os.W_OK):                         # running as root: permissions are not enforced
        pytest.skip("directory permissions not enforced")
    assert not busco_data.writable(readonly)
    assert busco_data.lineage_dir(str(readonly), local, "brassicales_odb12") == str(readonly)   # present: used
    assert busco_data.lineage_dir(str(readonly), local, "embryophyta_odb12") == local           # missing: results
    assert busco_data.compleasm_dir(str(readonly), local) == local                             # needs writing
    w = tmp_path / "writable"
    assert busco_data.lineage_dir(str(w), local, "embryophyta_odb12") == str(w)
    assert busco_data.compleasm_dir(str(w), local) == str(w)


def test_check_names(tmp_path):
    d = tmp_path / "db"
    assert busco_data.check_names(["eudicots_odb12"], [str(d)]) == []     # no local list: checked when downloading
    os.makedirs(d / "lineages")
    (d / "lineages" / "file_versions.tsv").write_text(VERSIONS)
    errs = busco_data.check_names(["eudicots_odb12", "brassicales_odb12"], [None, str(d)])
    assert len(errs) == 1 and "eudicots_odb12" in errs[0] and "eudicotyledons_odb12" in errs[0]
    dataset(d, "my_custom_odb12")                                          # present locally: accepted
    assert busco_data.check_names(["my_custom_odb12"], [str(d)]) == []


def test_fetch_present_needs_no_network(tmp_path, monkeypatch):
    dataset(tmp_path, "brassicales_odb12")
    monkeypatch.setattr(fetch_busco, "download", lambda *a, **k: pytest.fail("network used"))
    fetch_busco.ensure(str(tmp_path), ["brassicales_odb12"])
    assert (tmp_path / "lineages" / "brassicales_odb12.done").exists()      # compleasm marker


def test_fetch_unknown_name(tmp_path, monkeypatch):
    def fake_versions(ldir):
        (tmp_path / "lineages" / "file_versions.tsv").write_text(VERSIONS)
        return fetch_busco.read_versions(os.path.join(ldir, "file_versions.tsv"))
    monkeypatch.setattr(fetch_busco, "fetch_versions", fake_versions)
    with pytest.raises(fetch_busco.UnknownLineage, match="eudicotyledons_odb12"):
        fetch_busco.ensure(str(tmp_path), ["eudicots_odb12"])


def test_placement_file_names(tmp_path):
    (tmp_path / "v.tsv").write_text(VERSIONS)
    names = fetch_busco.placement_files(fetch_busco.read_versions(str(tmp_path / "v.tsv")))
    assert [(a, e) for a, e, _ in names] == [
        ("list_of_reference_markers.eukaryota_odb12.2025-01-15.txt.tar.gz",
         "list_of_reference_markers.eukaryota_odb12.2025-01-15.txt"),
        ("tree.eukaryota_odb12.2.2026-05-13.nwk.tar.gz", "tree.eukaryota_odb12.2.2026-05-13.nwk")]


def test_fetch_syncs_only_after_installing(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(fetch_busco.os, "sync", lambda: calls.append("sync"))
    dataset(tmp_path, "brassicales_odb12")
    fetch_busco.ensure(str(tmp_path), ["brassicales_odb12"])
    assert calls == []                                                 # nothing downloaded: no sync
    monkeypatch.setattr(fetch_busco, "fetch_versions", lambda ldir: {})
    monkeypatch.setattr(fetch_busco, "fetch_lineage", lambda ldir, name, rows: dataset(tmp_path, name))
    fetch_busco.ensure(str(tmp_path), ["embryophyta_odb12"])
    assert calls == ["sync"] and (tmp_path / "lineages" / "embryophyta_odb12.done").exists()
