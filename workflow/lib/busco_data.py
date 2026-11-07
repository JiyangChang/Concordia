"""Where the BUSCO datasets are read from / downloaded to (decided when the workflow starts).

- qc.busco_download_path set: a lineage already there is used from there; a missing one is downloaded there when the
  directory is writable, otherwise into the result directory (<outdir>/00_Data/busco_downloads).
- not set: everything goes to the result directory, shared by all BUSCO steps and BRAKER4.
- BRAKER4's compleasm refreshes file_versions.tsv.hash on every run, so its directory must be writable: the shared
  directory when writable, else the result directory.
"""
import difflib
import os


def writable(path):
    """True when `path` (or, if it does not exist yet, its nearest existing parent) can be written to."""
    p = os.path.abspath(path)
    while not os.path.exists(p):
        parent = os.path.dirname(p)
        if parent == p:
            return False
        p = parent
    return os.path.isdir(p) and os.access(p, os.W_OK | os.X_OK)


def has_lineage(directory, lineage):
    return bool(directory) and os.path.isfile(os.path.join(directory, "lineages", lineage, "dataset.cfg"))


def lineage_dir(shared, local, lineage):
    if shared and (has_lineage(shared, lineage) or writable(shared)):
        return shared
    return local


def compleasm_dir(shared, local):
    return shared if shared and writable(shared) else local


def check_names(lineages, directories):
    """Errors for lineage names missing from a local file_versions.tsv (none found locally: no check here; the
    download step checks against the online list)."""
    for d in directories:
        tsv = os.path.join(d, "lineages", "file_versions.tsv") if d else ""
        if tsv and os.path.isfile(tsv):
            break
    else:
        return []
    names = set()
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5 and f[4] == "lineages":
            names.add(f[0])
    errs = []
    for lin in lineages:
        if lin not in names and not any(has_lineage(d, lin) for d in directories):
            close = difflib.get_close_matches(lin, sorted(names), n=5, cutoff=0.6)
            errs.append(f"BUSCO dataset '{lin}' is not in {tsv}" +
                        (f" (similar names: {', '.join(close)})" if close else ""))
    return errs
