#!/usr/bin/env python3
"""GeMoMa-specific statistics: per-species contribution and evidence attributes.

GeMoMa runs once per homologous species (i=<species>, so every prediction ID starts with
"<species>_"); the per-species predictions are combined with GAF + AnnotationFinalizer.

--reference name=unfiltered.gff,species_final.gff (repeatable, one per species)
  unfiltered_predictions      GeMoMa predictions before filtering (o=true output)
  species_final_transcripts   transcripts in that species' own GeMoMa annotation
  merged_representative       merged transcripts whose representative prediction came from it
                              (its ID starts with "<species>_")
  merged_supported            merged transcripts predicted by this species at all: representative or
                              listed in the GAF 'alternative' attribute (identical CDS from another species)
--final: attributes of mRNA records in the merged set:
  evidence = number of species that predict this transcript,
  tie      = fraction of introns supported by RNA-seq (only with RNA evidence),
  tpc      = fraction of the transcript covered by RNA-seq,
  score    = GeMoMa alignment score.
"""
import argparse
import os
import statistics
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_gff, parse_attrs  # noqa: E402


def count_mrna(path):
    if not os.path.exists(path):
        return "-"
    return sum(1 for c in read_gff(path) if c[2] in ("mRNA", "prediction"))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def ids_of(d):
    """Prediction IDs that a merged transcript stands for: its own and its GAF alternatives."""
    own = [d.get("ID", ""), d.get("Name", ""), d.get("ref-gene", "")]
    alt = [x.strip('"') for x in d.get("alternative", "").split(",") if x]
    return own, alt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", required=True, help="merged GeMoMa annotation")
    ap.add_argument("--reference", action="append", required=True,
                    help="name=unfiltered_predictions.gff,species_final.gff")
    ap.add_argument("--out-references", required=True)
    ap.add_argument("--out-attributes", required=True)
    a = ap.parse_args()

    attrs = [parse_attrs(c[8]) for c in read_gff(a.final) if c[2] == "mRNA"]
    with open(a.out_references, "w") as fh:
        fh.write("species\tunfiltered_predictions\tspecies_final_transcripts\tmerged_representative\t"
                 "merged_supported\n")
        for item in a.reference:
            name, files = item.split("=", 1)
            unfiltered, species_final = (files.split(",", 1) + [""])[:2]
            pre = name + "_"
            rep = sup = 0
            for d in attrs:
                own, alt = ids_of(d)
                is_rep = any(x.startswith(pre) for x in own)
                rep += is_rep
                sup += is_rep or any(x.startswith(pre) for x in alt)
            fh.write(f"{name}\t{count_mrna(unfiltered)}\t{count_mrna(species_final) if species_final else '-'}\t"
                     f"{rep}\t{sup}\n")

    with open(a.out_attributes, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"merged_transcripts\t{len(attrs)}\n")
        ev = Counter(d.get("evidence") for d in attrs if d.get("evidence"))
        for k in sorted(ev, key=lambda x: num(x) or 0):
            fh.write(f"transcripts_supported_by_{k}_species\t{ev[k]}\n")
        for key, label in (("tie", "tie (intron RNA support)"), ("tpc", "tpc (RNA coverage)"), ("score", "score")):
            vals = [num(d.get(key)) for d in attrs]
            vals = [v for v in vals if v is not None and v == v]
            if not vals:
                fh.write(f"{label}\tnot available\n")
                continue
            fh.write(f"{label} mean\t{statistics.mean(vals):.3f}\n")
            if key in ("tie", "tpc"):
                fh.write(f"{label} = 1 (fully supported) %\t{100 * sum(v >= 1 for v in vals) / len(vals):.1f}\n")
                fh.write(f"{label} = 0 (no support) %\t{100 * sum(v <= 0 for v in vals) / len(vals):.1f}\n")


if __name__ == "__main__":
    main()
