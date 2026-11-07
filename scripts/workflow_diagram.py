#!/usr/bin/env python3
"""Turn a Snakemake rule graph (snakemake --rulegraph) into the documentation diagrams.

  workflow_diagram.py --rulegraph rg.dot --rules-dir workflow/rules --mode main --title "..." --out main.dot
  workflow_diagram.py --rulegraph rg.dot --rules-dir workflow/rules --mode full --out full.dot

main: only the main analysis steps, grouped by result directory (00_Data ... 06_Function) with a single
      07_Report node; QC / statistics steps are left out, and edges that ran through them are joined up
      (A -> qc -> B becomes A -> B). Every arrow is a direct dependency ("reads the output of"); arrows that
      are implied by a longer path are kept on purpose, so that e.g. all inputs of refine are visible.
full: every rule, also grouped by directory (QC / statistics steps in a light style).

The directory of a rule comes from its script("<dir>", ...) parameter in workflow/rules/*.smk; rules without
one fall back to the directory of their rule file.
"""
import argparse
import glob
import os
import re
from collections import defaultdict

DIRS = [("00_Data", "00_Data 输入与预处理 / Input", "#e8f1fb"),
        ("01_Repeat", "01_Repeat 重复序列 / Repeats", "#fdf0e3"),
        ("02_RNA", "02_RNA 比对与 EviAnn / RNA", "#e9f6ec"),
        ("03_Ab_initio", "03_Ab_initio BRAKER4", "#f3ecfa"),
        ("04_Homology", "04_Homology 同源注释 / Homology", "#fbeef0"),
        ("05_EVM", "05_EVM 整合·调优·修正·最终 / Integration", "#fff8dc"),
        ("06_Function", "06_Function 功能注释 / Function", "#eaf4f4"),
        ("07_Report", "07_Report 报告 / Report", "#f0f0f0")]
FILE_DIR = {"data": "00_Data", "assembly_qc": "00_Data", "repeat": "01_Repeat", "rna": "02_RNA",
            "abinitio": "03_Ab_initio", "homology": "04_Homology", "evm": "05_EVM", "qc": "05_EVM",
            "function": "06_Function", "report": "07_Report"}
# QC / statistics / bookkeeping rules: shown in the full graph only
NOT_MAIN = {"all", "report", "data_qc", "qc_summary", "assembly_stats", "assembly_qc_summary", "busco_genome",
            "tidk", "lai", "merqury", "fcs_gx", "nanoplot", "rna_align_summary", "reference_compare",
            "agat_config", "bam_sample_stats"}


def is_main(rule):
    return not (rule in NOT_MAIN or rule.startswith("qc_") or rule.endswith("_stats"))


def rule_dirs(rules_dir):
    """rule name (or name pattern for f-string rules) -> result directory."""
    out, patterns = {}, []
    for path in glob.glob(os.path.join(rules_dir, "*.smk")):
        default = FILE_DIR.get(os.path.basename(path)[:-4], "05_EVM")
        name = None
        for line in open(path):
            m = re.match(r"^\s*(?:local)?rule\s+(\w+)\s*:", line) or re.match(r'^\s*name:\s*f?"([^"]+)"', line)
            if m:
                name = m.group(1)
                if "{" in name:
                    patterns.append((re.compile("^" + re.sub(r"\\\{[^}]*\\\}", ".+", re.escape(name)) + "$"), default))
                else:
                    out[name] = default
                continue
            m = re.search(r'script\(\s*"(\d\d_\w+)"', line) or re.search(r"script\(\s*_step", line)
            if m and name and "{" not in name and m.groups():
                out[name] = m.group(1)
    return out, patterns


def read_rulegraph(path):
    labels, edges = {}, set()
    for line in open(path):
        m = re.match(r'\s*(\d+)\[label = "([^"]+)"', line)
        if m:
            labels[m.group(1)] = m.group(2)
        m = re.match(r"\s*(\d+) -> (\d+)", line)
        if m:
            edges.add((m.group(1), m.group(2)))
    return labels, {(labels[a], labels[b]) for a, b in edges}


def contract(nodes, edges, keep):
    """Remove nodes not in keep; connect their predecessors to their successors (transitively)."""
    succ = defaultdict(set)
    for a, b in edges:
        succ[a].add(b)
    out = set()
    for a in nodes:
        if a not in keep:
            continue
        stack, seen = list(succ[a]), set()
        while stack:
            b = stack.pop()
            if b in seen:
                continue
            seen.add(b)
            if b in keep:
                out.add((a, b))
            else:
                stack.extend(succ[b])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rulegraph", required=True)
    ap.add_argument("--rules-dir", required=True)
    ap.add_argument("--mode", choices=["main", "full"], required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    labels, edges = read_rulegraph(a.rulegraph)
    nodes = set(labels.values())
    named, patterns = rule_dirs(a.rules_dir)

    def where(r):
        if r in named:
            return named[r]
        return next((d for p, d in patterns if p.match(r)), "05_EVM")

    if a.mode == "main":
        keep = {n for n in nodes if is_main(n)}
        # every direct dependency is kept (no transitive reduction): an arrow means "reads the output of"
        edges = contract(nodes, edges, keep)
        nodes = keep
    groups = defaultdict(list)
    for n in sorted(nodes):
        if a.mode == "main" or n not in ("all",):
            groups[where(n) if n != "all" else "07_Report"].append(n)

    L = ["digraph workflow {",
         '  graph [rankdir=TB, fontname="Noto Sans CJK SC", fontsize=14, nodesep=0.25, ranksep=0.4, '
         f'labelloc=t, label="{a.title}", compound=true];',
         '  node [shape=box, style="rounded,filled", fillcolor=white, fontname="DejaVu Sans", fontsize=11];',
         '  edge [color="#8a8a8a", arrowsize=0.6];']
    for i, (d, title, colour) in enumerate(DIRS):
        members = groups.get(d, [])
        if a.mode == "main" and d == "07_Report":
            continue
        if not members:
            continue
        L.append(f'  subgraph cluster_{i} {{ label="{title}"; style="rounded,filled"; fillcolor="{colour}"; '
                 'color="#bbbbbb"; fontname="Noto Sans CJK SC";')
        for n in members:
            extra = "" if is_main(n) else ', fillcolor="#f7f7f7", fontcolor="#777777", color="#bbbbbb"'
            L.append(f'    "{n}" [label="{n}"{extra}];')
        L.append("  }")
    if a.mode == "main":
        L.append('  report [label="07_Report 图文报告 / report\\n（汇总各步骤统计与 BUSCO / OMArk / AGAT 质控）", '
                 'fillcolor="#f0f0f0", fontname="Noto Sans CJK SC"];')
        sinks = {n for n in nodes if not any(x == n for x, _ in edges)}
        for n in sorted(sinks):
            L.append(f'  "{n}" -> report [style=dashed];')
    for x, y in sorted(edges):
        L.append(f'  "{x}" -> "{y}";')
    L.append("}")
    with open(a.out, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"{a.mode}: {len(nodes)} rules, {len(edges)} edges -> {a.out}")


if __name__ == "__main__":
    main()
