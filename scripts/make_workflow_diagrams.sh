#!/usr/bin/env bash
# Regenerate the workflow diagrams in docs/ from the current Snakemake rules.
#   bash scripts/make_workflow_diagrams.sh            (run from the project directory)
# Needs snakemake and graphviz (dot) on PATH, e.g. in the "concordia" environment:
#   conda install -n concordia -c conda-forge graphviz
#
# The rule graph is built from the mock configuration tests/config.test.yaml (short reads + PacBio/ONT long
# reads, two homologous species, no external databases needed), with eviann.swissprot set as recommended.
# Outputs (docs/):
#   workflow_main.{dot,svg,png}        main steps grouped by result directory (QC / statistics left out)
#   workflow_main_norna.{dot,svg,png}  the same without RNA-seq (rna.allow_none: true)
#   workflow_rulegraph.{dot,svg,png}   every rule
set -euo pipefail
cd "$(dirname "$0")/.."
command -v snakemake >/dev/null || { echo "snakemake not found on PATH" >&2; exit 1; }
command -v dot >/dev/null || { echo "graphviz (dot) not found on PATH" >&2; exit 1; }

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
CFG=tests/config.test.yaml
SWISSPROT='eviann={swissprot: tests/data/orthodb.fa}'      # mock stand-in for UniProt-SwissProt
NORNA='rna={short: null, long: null, allow_none: true}'

rulegraph() {   # $1 = output dot, remaining = extra --config items
    local out=$1; shift
    snakemake -s workflow/Snakefile --configfile "$CFG" --nolock --quiet all --rulegraph --config "$@" > "$out"
}
render() {      # $1 = dot file without extension
    dot -Tsvg "$1.dot" -o "$1.svg"
    dot -Tpng -Gdpi=110 "$1.dot" -o "$1.png"
}

rulegraph "$TMP/rg.dot" "$SWISSPROT"
rulegraph "$TMP/rg_norna.dot" "$NORNA"

python3 scripts/workflow_diagram.py --rulegraph "$TMP/rg.dot" --rules-dir workflow/rules --mode main \
    --title "主流程（有 RNA-seq：短读 + 长读）/ Main workflow (with RNA-seq)" --out docs/workflow_main.dot
python3 scripts/workflow_diagram.py --rulegraph "$TMP/rg_norna.dot" --rules-dir workflow/rules --mode main \
    --title "主流程（无 RNA 模式，rna.allow_none: true）/ Main workflow without RNA-seq" --out docs/workflow_main_norna.dot
python3 scripts/workflow_diagram.py --rulegraph "$TMP/rg.dot" --rules-dir workflow/rules --mode full \
    --title "全部步骤（有 RNA-seq）/ All rules (with RNA-seq)" --out docs/workflow_rulegraph.dot

for f in docs/workflow_main docs/workflow_main_norna docs/workflow_rulegraph; do
    render "$f"
    echo "written $f.{dot,svg,png}"
done
