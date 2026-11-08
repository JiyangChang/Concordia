#!/usr/bin/env python3
"""Make sure BUSCO datasets are present in a download directory; fetch the missing ones.

  fetch_busco.py --dir <busco_downloads> --lineage brassicales_odb12 [--lineage ...] [--compleasm] [--done flag]

Layout and files are the same as BUSCO's / compleasm's own downloaders and scripts/download_databases.sh:
  <dir>/lineages/<lineage>/           extracted dataset (dataset.cfg marks a complete one)
  <dir>/lineages/<lineage>.done       marker compleasm checks (without it compleasm downloads the lineage again)
  --compleasm additionally: file_versions.tsv(.hash, .done), eukaryota_odb12 and the eukaryota placement files,
  which BRAKER4's compleasm needs.
Datasets already present are not touched (no network access is needed then). Downloads are checked against the md5
in file_versions.tsv and moved into place only when complete, so parallel runs sharing the directory are safe.
"""
import argparse
import difflib
import hashlib
import os
import shutil
import sys
import tarfile
import time
import urllib.request

BASE = "https://busco-data.ezlab.org/v5/data"
COMPLEASM_LINEAGE = "eukaryota_odb12"


class FetchError(Exception):
    pass


class UnknownLineage(FetchError):
    pass


def log(msg):
    print(f"[{time.strftime('%F %T')}] {msg}", flush=True)


def tmp_name(path):
    return f"{path}.part{os.getpid()}"


def download(url, dest, tries=3):
    tmp = tmp_name(dest)
    for i in range(1, tries + 1):
        try:
            log(f"downloading {url}")
            with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as out:
                shutil.copyfileobj(r, out, 1 << 20)
            os.replace(tmp, dest)
            return
        except OSError as e:          # URLError, timeouts, connection resets
            if os.path.exists(tmp):
                os.remove(tmp)
            if i == tries:
                raise FetchError(f"cannot download {url}: {e}")
            log(f"  failed ({e}); retrying in 20 s")
            time.sleep(20)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_versions(path):
    """name -> (date, md5, type) from file_versions.tsv."""
    rows = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5:
            rows[(f[0], f[4])] = (f[1], f[2])
    return rows


def suggest(name, rows):
    names = sorted({n for n, t in rows if t == "lineages"})
    close = difflib.get_close_matches(name, names, n=5, cutoff=0.6)
    return f" (similar names: {', '.join(close)})" if close else ""


def fetch_versions(ldir):
    """Current file_versions.tsv (+ .hash) into <dir>/lineages; md5-checked."""
    tsv, hsh = os.path.join(ldir, "file_versions.tsv"), os.path.join(ldir, "file_versions.tsv.hash")
    tmp_tsv, tmp_hsh = tmp_name(tsv) + ".tsv", tmp_name(hsh) + ".hash"
    download(f"{BASE}/file_versions.tsv.hash", tmp_hsh)
    download(f"{BASE}/file_versions.tsv", tmp_tsv)
    want = open(tmp_hsh).read().split()[0]
    if md5(tmp_tsv) != want:
        os.remove(tmp_tsv)
        os.remove(tmp_hsh)
        raise FetchError("md5 of file_versions.tsv does not match file_versions.tsv.hash")
    os.replace(tmp_tsv, tsv)
    os.replace(tmp_hsh, hsh)
    return read_versions(tsv)


def fetch_archive(url, want_md5, ldir):
    """Download a .tar.gz, check its md5 and extract it into a temporary directory inside <ldir>; returns that dir."""
    arc = os.path.join(ldir, tmp_name(os.path.basename(url)))
    download(url, arc)
    try:
        got = md5(arc)
        if got != want_md5:
            raise FetchError(f"md5 of {os.path.basename(url)} is {got}, expected {want_md5}")
        tmp = os.path.join(ldir, f".extract{os.getpid()}")
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
        with tarfile.open(arc) as t:
            t.extractall(tmp, filter="data") if hasattr(tarfile, "data_filter") else t.extractall(tmp)
        return tmp
    finally:
        if os.path.exists(arc):
            os.remove(arc)


def install(tmp, name, dest_dir):
    """Move <tmp>/<name> to <dest_dir>/<name> unless another run has put it there meanwhile."""
    src, dst = os.path.join(tmp, name), os.path.join(dest_dir, name)
    if not os.path.exists(src):
        raise FetchError(f"archive does not contain {name}")
    if not os.path.exists(dst):
        os.replace(src, dst)
    shutil.rmtree(tmp, ignore_errors=True)


def has_lineage(ldir, name):
    return os.path.isfile(os.path.join(ldir, name, "dataset.cfg"))


def fetch_lineage(ldir, name, rows):
    if (name, "lineages") not in rows:
        raise UnknownLineage(f"BUSCO dataset '{name}' does not exist in {BASE}/file_versions.tsv{suggest(name, rows)}")
    date, want = rows[(name, "lineages")]
    tmp = fetch_archive(f"{BASE}/lineages/{name}.{date}.tar.gz", want, ldir)
    install(tmp, name, ldir)
    log(f"{name} ({date}) installed")


def placement_files(rows):
    """(archive name, extracted file name, md5) of every eukaryota placement file, as compleasm downloads them."""
    out = []
    for (name, typ), (date, want) in sorted(rows.items()):
        if typ == "placement_files" and "eukaryota" in name:
            p = name.split(".")
            extracted = ".".join(p[:-1] + [date, p[-1]])
            out.append((extracted + ".tar.gz", extracted, want))
    return out


def ensure(dirpath, lineages, compleasm=False):
    ldir = os.path.join(dirpath, "lineages")
    os.makedirs(ldir, exist_ok=True)
    wanted = list(dict.fromkeys(lineages + ([COMPLEASM_LINEAGE] if compleasm else [])))
    missing = [l for l in wanted if not has_lineage(ldir, l)]
    pdir = os.path.join(ldir, "placement_files")
    need_place = compleasm and not os.path.exists(os.path.join(ldir, "placement_files.done"))
    for l in wanted:
        if l not in missing:
            log(f"{l}: present in {ldir}")
    rows = None
    if missing or need_place or (compleasm and not os.path.exists(os.path.join(ldir, "file_versions.tsv"))):
        rows = fetch_versions(ldir)
        for l in missing:
            fetch_lineage(ldir, l, rows)
    if need_place:
        os.makedirs(pdir, exist_ok=True)
        for arc, extracted, want in placement_files(rows):
            if os.path.exists(os.path.join(pdir, extracted)):
                continue
            tmp = fetch_archive(f"{BASE}/placement_files/{arc}", want, pdir)
            install(tmp, extracted, pdir)
        open(os.path.join(ldir, "placement_files.done"), "w").close()
        log("eukaryota placement files installed")
    if compleasm:
        open(os.path.join(ldir, "file_versions.tsv.done"), "a").close()
    if missing or need_place:
        # flush to disk now: on removable / network storage a disconnect could otherwise lose data that was
        # already reported as installed (seen with a USB disk: the dataset vanished after re-mounting)
        os.sync()
    for l in wanted:                                   # compleasm's per-lineage markers
        marker = os.path.join(ldir, f"{l}.done")
        if not os.path.exists(marker):
            try:
                open(marker, "w").close()
            except OSError:                            # read-only shared directory with the dataset already there
                pass


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="BUSCO download directory (holds lineages/)")
    ap.add_argument("--lineage", action="append", default=[])
    ap.add_argument("--compleasm", action="store_true", help="also the files BRAKER4's compleasm needs")
    ap.add_argument("--done", help="flag file written when everything is present")
    a = ap.parse_args()
    d = os.path.abspath(a.dir)
    try:
        ensure(d, a.lineage, a.compleasm)
    except UnknownLineage as e:
        sys.stderr.write(f"ERROR: {e}\nCheck qc.busco_lineages.\n")
        sys.exit(1)
    except FetchError as e:
        cmd = " ".join([f"python3 {os.path.abspath(__file__)} --dir {d}"] + [f"--lineage {l}" for l in a.lineage] +
                       (["--compleasm"] if a.compleasm else []))
        sys.stderr.write(f"ERROR: {e}\nIf this node has no internet access, run on a machine that has (same paths):\n"
                         f"  {cmd}\nor download the dataset as described in docs/deployment_guide.md, then re-run.\n")
        sys.exit(1)
    if a.done:
        os.makedirs(os.path.dirname(os.path.abspath(a.done)), exist_ok=True)
        with open(a.done, "w") as f:
            f.write(f"{d}\t{','.join(a.lineage)}{chr(9) + 'compleasm' if a.compleasm else ''}\n")


if __name__ == "__main__":
    main()
