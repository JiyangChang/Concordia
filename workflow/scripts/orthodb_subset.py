#!/usr/bin/env python3
"""Extract the proteins of selected OrthoDB taxonomic levels from an OrthoDB protein FASTA.

The BRAKER OrthoDB partitions (e.g. Viridiplantae.fa.gz) name each protein <species_id>:<gene>,
e.g. ">1000413_0:000000". OrthoDB's level2species table maps every species_id to its taxonomic
path, so a level (by name, e.g. Brassicales, or NCBI taxid, e.g. 3699) selects all species below it.

  orthodb_subset.py --fasta Viridiplantae.fa.gz --levels odb12v2_levels.tab.gz \
      --level2species odb12v2_level2species.tab.gz --select Brassicales Malvales \
      --out subset.faa --report subset_report.tsv

--list prints the levels available in the FASTA (with species / protein counts) and exits.
"""
import argparse
import gzip
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def to_int(x):
    try:
        return int(x)
    except ValueError:
        return 0


def read_levels(path):
    """taxid -> (name, genes, species)"""
    out = {}
    with opener(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 5:
                out[t[0]] = (t[1], to_int(t[2]), to_int(t[4]))
    return out


def read_level2species(path):
    """species_id -> list of level taxids (root ... species)"""
    out = {}
    with opener(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 4:
                out[t[1]] = t[3].strip("{}").split(",")
    return out


def species_in_fasta(path):
    seen = {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                sp = line[1:].split(":", 1)[0]
                seen[sp] = seen.get(sp, 0) + 1
    return seen


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="OrthoDB protein FASTA (.gz allowed)")
    ap.add_argument("--levels", required=True, help="odb12v*_levels.tab(.gz)")
    ap.add_argument("--level2species", required=True, help="odb12v*_level2species.tab(.gz)")
    ap.add_argument("--select", nargs="+", help="level names or NCBI taxids (union)")
    ap.add_argument("--list", action="store_true", help="list the levels present in --fasta and exit")
    ap.add_argument("--out", help="output protein FASTA")
    ap.add_argument("--report", help="report TSV")
    a = ap.parse_args()

    levels = read_levels(a.levels)
    l2s = read_level2species(a.level2species)
    by_name = {v[0].lower(): k for k, v in levels.items()}

    if a.list:
        counts = species_in_fasta(a.fasta)
        per = {}
        for sp, n in counts.items():
            for tax in l2s.get(sp, []):
                s = per.setdefault(tax, [0, 0])
                s[0] += 1
                s[1] += n
        print("taxid\tlevel\tspecies_in_fasta\tproteins_in_fasta")
        for tax, (nsp, nprot) in sorted(per.items(), key=lambda kv: -kv[1][1]):
            if tax in levels and nsp >= 2:
                print(f"{tax}\t{levels[tax][0]}\t{nsp}\t{nprot}")
        unmapped = [sp for sp in counts if sp not in l2s]
        if unmapped:
            print(f"# {len(unmapped)} species in the FASTA are not in level2species (e.g. {unmapped[0]})",
                  file=sys.stderr)
        return

    if not (a.select and a.out):
        ap.error("--select and --out are required unless --list is given")
    wanted, names = set(), []
    for s in a.select:
        tax = s if s in levels else by_name.get(str(s).lower())
        if tax is None:
            sys.exit(f"ERROR: OrthoDB level '{s}' not found in {a.levels} (use --list to see the choices)")
        wanted.add(tax)
        names.append(f"{levels[tax][0]} ({tax})")
    keep_species = {sp for sp, path in l2s.items() if wanted & set(path)}

    kept_sp, n_in, n_out, keep = set(), 0, 0, False
    with opener(a.fasta) as fh, open(a.out, "w") as out:
        for line in fh:
            if line.startswith(">"):
                n_in += 1
                sp = line[1:].split(":", 1)[0]
                keep = sp in keep_species
                if keep:
                    n_out += 1
                    kept_sp.add(sp)
            if keep:
                out.write(line)
    if n_out == 0:
        sys.exit(f"ERROR: no protein of {', '.join(names)} found in {a.fasta}")
    report = (f"selected_levels\t{'; '.join(names)}\n"
              f"species_in_levels\t{len(keep_species)}\n"
              f"species_in_subset\t{len(kept_sp)}\n"
              f"proteins_in_fasta\t{n_in}\n"
              f"proteins_in_subset\t{n_out}\n")
    if a.report:
        with open(a.report, "w") as fh:
            fh.write("item\tvalue\n" + report)
    sys.stdout.write(report)


if __name__ == "__main__":
    main()
