#!/usr/bin/env python3
"""Rename genome sequence IDs only when needed (EDTA: IDs <= 13 chars, [A-Za-z0-9_] only).

If ANY ID is longer than --max-len (or contains invalid characters and --on-invalid-chars is
set), ALL sequences are renamed to <prefix><zero-padded index> in input order, so new names can
never collide with original ones. Otherwise names are kept. A seqid_map.tsv (new<TAB>original)
is always written so downstream steps can treat both cases the same way.
Sequences shorter than --min-len are dropped (recorded in the report).
With --link-unchanged, an input that needs no change (no renaming, nothing dropped) is linked
instead of copied (saves disk space for large genomes).
"""
import argparse
import re
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_fasta, write_fasta  # noqa: E402

VALID = re.compile(r"^[A-Za-z0-9_]+$")


def decide(ids, max_len, on_invalid):
    too_long = [i for i in ids if len(i) > max_len]
    invalid = [i for i in ids if not VALID.match(i)]
    if too_long:
        return True, f"{len(too_long)} ID(s) longer than {max_len} chars (e.g. {too_long[0]})", invalid
    if invalid:
        if on_invalid:
            return True, f"{len(invalid)} ID(s) with characters outside [A-Za-z0-9_] (e.g. {invalid[0]})", invalid
        raise SystemExit(f"ERROR: IDs with invalid characters for EDTA (e.g. {invalid[0]}) and "
                         "rename_seqids.on_invalid_chars is false")
    return False, f"all IDs <= {max_len} chars and valid; names kept", invalid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--max-len", type=int, default=13)
    ap.add_argument("--prefix", default="seq")
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--on-invalid-chars", choices=["true", "false"], default="true")
    ap.add_argument("--link-unchanged", action="store_true",
                    help="symlink --out-fasta to the input when no sequence is renamed or dropped")
    a = ap.parse_args()

    records = []
    dropped = 0
    for header, seq in read_fasta(a.fasta):
        sid = header.split()[0]
        if len(seq) < a.min_len:
            dropped += 1
            continue
        records.append((sid, seq))
    if not records:
        raise SystemExit(f"ERROR: no sequences left after removing {dropped} sequence(s) shorter than "
                         f"--min-len {a.min_len}")
    ids = [r[0] for r in records]
    if len(set(ids)) != len(ids):
        raise SystemExit("ERROR: duplicate sequence IDs in genome FASTA")

    rename, reason, _ = decide(ids, a.max_len, a.on_invalid_chars == "true")
    width = max(4, len(str(len(records))))
    if rename and len(a.prefix) + width > a.max_len:
        raise SystemExit(f"ERROR: prefix '{a.prefix}' + {width} digits exceeds {a.max_len} chars")

    link = a.link_unchanged and not rename and not dropped
    if os.path.lexists(a.out_fasta):
        os.remove(a.out_fasta)
    if link:
        os.symlink(os.path.realpath(a.fasta), a.out_fasta)
    with open(a.out_map, "w") as mp:
        mp.write("new_id\toriginal_id\n")
        fa = None if link else open(a.out_fasta, "w")
        for i, (sid, seq) in enumerate(records, 1):
            new = f"{a.prefix}{i:0{width}d}" if rename else sid
            if fa:
                write_fasta(fa, new, seq)
            mp.write(f"{new}\t{sid}\n")
        if fa:
            fa.close()
    with open(a.out_report, "w") as rp:
        rp.write(f"renamed\t{'yes' if rename else 'no'}\n")
        rp.write(f"reason\t{reason}\n")
        rp.write(f"sequences_kept\t{len(records)}\n")
        rp.write(f"sequences_dropped_min_len\t{dropped}\n")
        rp.write(f"fasta\t{'symbolic link to the input (unchanged)' if link else 'written'}\n")


if __name__ == "__main__":
    main()
