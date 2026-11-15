#!/usr/bin/env python3
"""07_Report: build a self-contained HTML report (text + tables + figures) from all steps, in Chinese and English.

  make_report.py --outdir <results> ... --out-zh annotation_report_zh.html --out-en annotation_report_en.html

Both versions come from the same data in one run; all text is given as L("中文", "English") pairs.
Each section is skipped gracefully when its input does not exist.
"""
import argparse
import textwrap
import base64
import csv
import glob
import html
import io
import os
import re
from collections import Counter

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # report still renders, without figures
    plt = None

LANG = "zh"     # language of the report being built (set per output in main)


def L(zh, en):
    """Text in the language of the report being built."""
    return zh if LANG == "zh" else en


CSS = """
:root{--bg:#fff;--fg:#1d2433;--muted:#5b6475;--line:#d9dde5;--accent:#2f6fb0;--card:#f6f8fb}
@media (prefers-color-scheme: dark){:root{--bg:#14171d;--fg:#e6e9ef;--muted:#9aa3b2;--line:#2c323d;--accent:#7fb0e6;--card:#1b1f27}}
body{background:var(--bg);color:var(--fg);font:15px/1.6 system-ui,-apple-system,"Noto Sans CJK SC",sans-serif;margin:0}
.layout{display:grid;grid-template-columns:250px minmax(0,1fr);max-width:1420px;margin:0 auto}
nav.side{position:sticky;top:0;align-self:start;height:100vh;overflow-y:auto;padding:24px 10px 24px 16px;
  border-right:1px solid var(--line);box-sizing:border-box}
nav.side .navtitle{font-weight:600;margin:0 0 10px 8px}
nav.side a{display:block;padding:5px 8px;margin:1px 0;border-radius:5px;color:var(--fg);text-decoration:none;
  font-size:13.5px;line-height:1.35;border-left:3px solid transparent}
nav.side a:hover{background:var(--card)}
nav.side a.active{background:var(--card);color:var(--accent);border-left-color:var(--accent);font-weight:600}
nav.side .step{display:block;font-size:11px;color:var(--muted);font-weight:400}
nav.side .top{margin-top:12px;font-size:12px;color:var(--muted)}
nav.side .lang{margin:0 0 12px 8px;font-size:12px}
main{max-width:1100px;margin:0 auto;padding:24px 24px 64px;min-width:0}
section{scroll-margin-top:12px}
@media (max-width:900px){.layout{display:block}
  nav.side{position:static;height:auto;border-right:0;border-bottom:1px solid var(--line);padding:12px 16px;
    display:flex;flex-wrap:wrap;gap:4px} nav.side .navtitle{width:100%;margin:0 0 4px}
  nav.side a{border-left:0;padding:3px 8px;background:var(--card)} nav.side .step{display:inline;margin-right:4px}
  main{padding:16px 16px 64px}}
h1{font-size:26px;margin:0 0 4px} h2{font-size:20px;border-bottom:1px solid var(--line);padding-bottom:4px;margin-top:40px}
.muted{color:var(--muted)}
.tbl{overflow-x:auto;margin:12px 0}
table{border-collapse:collapse;margin:0 auto;font-size:13px}
th,td{border:1px solid var(--line);padding:4px 8px;text-align:left;white-space:nowrap} th{background:var(--card)}
img{max-width:100%;height:auto;background:#fff;border-radius:6px;margin:8px auto;display:block}
details{margin:8px 0} summary{cursor:pointer} details table{width:100%} details td{white-space:normal} details td:first-child,details td:nth-child(2){white-space:nowrap} h4{margin:12px 0 4px}
.note{background:var(--card);border-left:3px solid var(--accent);padding:8px 12px;margin:12px 0}
figure{margin:12px auto 20px;text-align:center}
figcaption{font-size:13px;color:var(--muted);margin:2px auto 0;line-height:1.5;max-width:860px;text-align:left}
figcaption b{color:var(--fg)}
"""
FUNC_SRC_COL = {"InterPro": "#2f6fb0", "Pfam": "#2f6fb0", "COG_category": "#e0a43a", "KEGG_ko": "#e0a43a",
                "KEGG_Pathway": "#e0a43a", "EC": "#e0a43a", "Preferred_name": "#e0a43a", "Description": "#e0a43a"}
SOURCE_LABEL = {"eviann": "EviAnn", "braker4": "BRAKER4", "gemoma": "GeMoMa_merged", "evm": "EVM",
                "refine": "refine", "final": "final"}


def src_label(key):
    return SOURCE_LABEL.get(key, key.replace("gemoma_", "GeMoMa_", 1))


# highlight the section currently in view in the side navigation (no external libraries)
NAV_JS = """
(function(){
  var links = {};
  document.querySelectorAll('nav.side a[href^="#"]').forEach(function(a){ links[a.getAttribute('href').slice(1)] = a; });
  var secs = Array.prototype.slice.call(document.querySelectorAll('main section'));
  function update(){
    var cur = secs[0];
    secs.forEach(function(s){ if (s.getBoundingClientRect().top <= 80) cur = s; });
    Object.keys(links).forEach(function(k){ links[k].classList.toggle('active', cur && k === cur.id); });
    var a = cur && links[cur.id];
    if (a && a.scrollIntoView && window.innerWidth > 900) {
      var r = a.getBoundingClientRect(), n = a.parentNode.getBoundingClientRect();
      if (r.top < n.top || r.bottom > n.bottom) a.scrollIntoView({block: 'nearest'});
    }
  }
  window.addEventListener('scroll', update, {passive: true});
  update();
})();
"""

FIG_MARK = "@@FIG@@"   # replaced by the figure number in document order (see number_figures)


def read_tsv(path):
    if not path or not os.path.exists(path):
        return None, []
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("##")]
    if not rows:
        return None, []
    return rows[0], rows[1:]


def read_tsv_ws(path):
    """Whitespace-separated table with a header line (e.g. LAI output)."""
    rows = [l.split() for l in open(path) if l.strip()]
    return (rows[0], rows[1:]) if rows else (None, [])


def table(header, rows, limit=200):
    if not header:
        return "<p class='muted'>" + L("（无数据）", "(no data)") + "</p>"
    h = "".join(f"<th>{html.escape(c.lstrip('#'))}</th>" for c in header)
    b = "".join("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in r) + "</tr>" for r in rows[:limit])
    more = ("<p class='muted'>" + L(f"仅显示前 {limit} 行，共 {len(rows)} 行。",
                                    f"Only the first {limit} of {len(rows)} rows are shown.") + "</p>"
            if len(rows) > limit else "")
    return f"<div class='tbl'><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>{more}"


def fig_to_img(fig, figdir, name, caption=""):
    """Embed a figure; caption = (title, explanation) or text, shown as "Figure N. title explanation".
    Figures carry English labels only, so both language versions share the PNG files in figures/."""
    os.makedirs(figdir, exist_ok=True)
    fig.savefig(os.path.join(figdir, name + ".png"), dpi=150, bbox_inches="tight")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    img = f"<img alt='{name}' src='data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}'>"
    title, text = caption if isinstance(caption, tuple) else ("", caption)
    lab = L("图", "Figure")
    paren = L(f"（{html.escape(name)}.png）", f"({html.escape(name)}.png)")
    cap = (f"<figcaption><b>{lab} {FIG_MARK}. {html.escape(title)}</b> {html.escape(text)} "
           f"<span class='muted'>{paren}</span></figcaption>")
    return f"<figure>{img}{cap}</figure>"


def legend_below(ax, handles=None, ncol=4, gap_in=0.62):
    """Legend centred below the x-axis label: offset in inches, independent of the figure height."""
    fig = ax.figure
    ax_h = ax.get_position().height * fig.get_figheight()
    kw = dict(ncol=ncol, loc="upper center", bbox_to_anchor=(0.5, -gap_in / ax_h), frameon=False)
    return ax.legend(handles=handles, **kw) if handles else ax.legend(**kw)


def number_figures(doc):
    """Number figure captions in the order they appear in the final page."""
    parts = doc.split(FIG_MARK)
    return "".join(p + (str(i + 1) if i < len(parts) - 1 else "") for i, p in enumerate(parts))


def busco_figure(header, rows, figdir, name="busco_comparison", caption=""):
    if plt is None or not header:
        return ""
    idx = {c: i for i, c in enumerate(header)}
    rows = [r for r in rows if r[idx["S"]]]
    if not rows:
        return ""
    labels = [f"{src_label(r[idx['source']])}\n{r[idx['lineage']]}" for r in rows]
    parts = [("S", "Single", "#2f6fb0"), ("D", "Duplicated", "#7fb0e6"), ("F", "Fragmented", "#e0a43a"), ("M", "Missing", "#c8553d")]
    fig, ax = plt.subplots(figsize=(8, 0.5 * len(rows) + 1.2))
    left = [0.0] * len(rows)
    for key, part_label, color in parts:      # not "name": that is the figure file name
        vals = [float(r[idx[key]]) for r in rows]
        ax.barh(labels, vals, left=left, color=color, label=part_label)
        left = [l + v for l, v in zip(left, vals)]
    for i, r in enumerate(rows):                       # complete (S + D) at the bar end
        c_ = float(r[idx["S"]]) + float(r[idx["D"]])
        ax.text(101, i, f"C {c_:.1f}%", va="center", fontsize=8)
    ax.set_xlim(0, 112); ax.set_xticks(range(0, 101, 20)); ax.set_xlabel("BUSCO (%)"); ax.invert_yaxis()
    legend_below(ax, ncol=4)
    return fig_to_img(fig, figdir, name, caption or (
        L("BUSCO 评估。", "BUSCO assessment."),
        L("每行一个注释来源（或基因组）× 一个 BUSCO 数据集；条形按单拷贝完整（S）、重复完整（D）、"
          "片段（F）、缺失（M）的百分比堆叠，右侧为完整度 C = S + D。",
          "One row per annotation source (or the genome) x BUSCO dataset; bars stack the percentages of complete "
          "single-copy (S), complete duplicated (D), fragmented (F) and missing (M) BUSCOs; completeness C = S + D "
          "is shown on the right.")))


def _omark_parts(text):
    """'A:88.07%[P:5.77%,F:1.44%],I:4.06%[...],...' -> {'A': 88.07, 'A.P': 5.77, 'A.F': 1.44, 'I': 4.06, ...}"""
    out = {}
    for m in re.finditer(r"([A-Z]):([\d.]+)%(?:\[([^\]]*)\])?", text or ""):
        out[m.group(1)] = float(m.group(2))
        for k, v in re.findall(r"([A-Z]):([\d.]+)%", m.group(3) or ""):
            out[f"{m.group(1)}.{k}"] = float(v)
    return out


def omark_figure(header, rows, figdir, name="omark_comparison"):
    """OMArk of every annotation source: completeness (S / D / M of the conserved HOGs) and consistency (A split into
    full-length / partial / fragment, I, C, U of all proteins), as stacked bars like the BUSCO figure."""
    if plt is None or not header or "omark_completeness" not in header:
        return ""
    ic, ik = header.index("omark_completeness"), header.index("omark_consistency")
    seen, data = set(), []
    for r in rows:                                  # one row per source (summary.tsv repeats it per lineage)
        if r[0] in seen or not r[ic].strip():
            continue
        seen.add(r[0])
        data.append((src_label(r[0]), _omark_parts(r[ic]), _omark_parts(r[ik])))
    if not data:
        return ""
    comp = [("S", "Single", "#2f6fb0"), ("D", "Duplicated", "#7fb0e6"), ("M", "Missing", "#c8553d")]
    cons = [("A_full", "Consistent", "#5a9e6f"), ("A.P", "Consistent partial", "#8cc49c"),
            ("A.F", "Consistent fragment", "#c4e3cc"), ("I", "Inconsistent", "#e0a43a"),
            ("C", "Contaminant", "#8c6bb1"), ("U", "Unknown", "#999999")]
    for _, _, c in data:
        c["A_full"] = max(0.0, c.get("A", 0) - c.get("A.P", 0) - c.get("A.F", 0))
    labels = [d[0] for d in data]
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.45 * len(data) + 1.6), sharey=True)
    for ax, parts, idx, title in ((axes[0], comp, 1, "Completeness (conserved HOGs)"),
                                  (axes[1], cons, 2, "Consistency (all proteins)")):
        left = [0.0] * len(data)
        for key, lab, col in parts:
            vals = [d[idx].get(key, 0.0) for d in data]
            ax.barh(labels, vals, left=left, color=col, label=lab)
            left = [l + v for l, v in zip(left, vals)]
        ax.set_xlim(0, 100); ax.set_xlabel("%"); ax.set_title(title, fontsize=10)
        legend_below(ax, ncol=3 if len(parts) > 3 else len(parts))
    axes[0].invert_yaxis()
    for i, d in enumerate(data):                    # values at the bar end of the main categories
        axes[0].text(101, i, f"S+D {d[1].get('S', 0) + d[1].get('D', 0):.1f}%", va="center", fontsize=8)
    axes[0].set_xlim(0, 118)
    return fig_to_img(fig, figdir, name, (
        L("各注释来源的 OMArk 评估。", "OMArk assessment of the annotation sources."),
        L("左：完整度，即该谱系保守的 HOG（直系同源组）中被找到的比例，分单拷贝（S）、多拷贝（D）、缺失（M）；"
          "右：一致性，即全部蛋白中与该物种所在谱系一致（A，又分全长、部分 partial、片段 fragment）、谱系不一致（I）、"
          "疑似污染（C）和找不到同源（U）的比例。各来源每个基因用最长的蛋白。",
          "Left: completeness, the share of the conserved HOGs (hierarchical orthologous groups) of the lineage that "
          "are found, as single-copy (S), duplicated (D) or missing (M); right: consistency, the share of all "
          "proteins that are consistent with the species' lineage (A, split into full-length, partial and fragment), "
          "inconsistent (I), likely contaminants (C) or without homologs (U). The longest protein of every gene is used.")))


def bar_figure(labels, values, xlabel, figdir, name, caption="", colors=None, legend=None, wrap=0, values_text=None):
    """Horizontal bars. colors: one colour per bar; legend: [(label, colour)]; wrap: wrap labels at n chars;
    values_text: format string for a value label at the bar end (e.g. "{:.1f}")."""
    if plt is None or not labels:
        return ""
    if wrap:
        labels = [textwrap.fill(l_, wrap) for l_ in labels]
    fig, ax = plt.subplots(figsize=(7.5, 0.4 * len(labels) + 1))
    ax.barh(range(len(labels)), values, color=colors or "#2f6fb0")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis(); ax.set_xlabel(xlabel)
    if values_text:
        top = max(values) if values else 0
        for i, v in enumerate(values):
            ax.text(v + top * 0.01, i, values_text.format(v), va="center", fontsize=8)
        ax.set_xlim(0, top * 1.12 if top else 1)
    if legend:
        from matplotlib.patches import Patch
        legend_below(ax, handles=[Patch(color=c, label=l_) for l_, c in legend], ncol=len(legend))
    return fig_to_img(fig, figdir, name, caption)


def kv_file(path):
    out = []
    if path and os.path.exists(path):
        for line in open(path):
            if "\t" in line:
                out.append(line.rstrip("\n").split("\t", 1))
    return out


def protein_qc_rows(qc_dir):
    rows = []
    for f_ in sorted(glob.glob(os.path.join(qc_dir, "busco", "*", "short_summary.txt"))):
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", open(f_).read())
        if m:
            rows.append(["BUSCO " + os.path.basename(os.path.dirname(f_)),
                         "C:{}% [S:{}%, D:{}%], F:{}%, M:{}%, n:{}".format(*m.groups())])
    f_ = os.path.join(qc_dir, "omark", "omark.sum")
    if os.path.exists(f_):
        for line in open(f_):
            if line.startswith(("S:", "A:")):
                rows.append(["OMArk " + ("completeness" if line.startswith("S:") else "consistency"), line.strip()])
    return rows


def distribution_figure(dist, figdir, name, what=""):
    """Gene-structure distributions of one or more annotation sets, as % of each set so that sets of
    different size can be compared: isoforms per gene (grouped bars 1, 2, 3, 4, >=5; omitted when every
    set has one transcript per gene), exons per transcript and CDS length (step outlines)."""
    if plt is None or not dist:
        return ""
    series = {}
    for st_, kind, v in dist:
        series.setdefault((st_, kind), []).append(float(v))
    sets = list(dict.fromkeys(k[0] for k in series))
    pal = ["#2f6fb0", "#e0a43a", "#5a9e6f", "#c8553d", "#8c6bb1", "#7fb0e6", "#999999", "#d4a5a5"]
    colors = {st_: pal[i % len(pal)] for i, st_ in enumerate(sets)}   # many homolog species -> cycle
    iso_sets = [st_ for st_ in sets if any(v > 1 for v in series.get((st_, "isoforms_per_gene"), []))]
    panels = (["iso"] if iso_sets else []) + ["exons", "cds"]
    f, axes = plt.subplots(1, len(panels), figsize=(3.8 * len(panels), 3.4))
    axes = list(axes) if len(panels) > 1 else [axes]
    for ax, panel in zip(axes, panels):
        if panel == "iso":
            cats = ["1", "2", "3", "4", "≥5"]
            w = 0.8 / len(sets)
            for i, st_ in enumerate(sets):
                vals = series.get((st_, "isoforms_per_gene"), [])
                if not vals:
                    continue
                cnt = [sum(1 for v in vals if int(v) == k) for k in (1, 2, 3, 4)] + [sum(1 for v in vals if v >= 5)]
                ax.bar([j + (i - (len(sets) - 1) / 2) * w for j in range(5)], [100 * c / len(vals) for c in cnt],
                       width=w, color=colors[st_], label=st_)
            ax.set_xticks(range(5)); ax.set_xticklabels(cats)
            ax.set_xlabel("transcripts per gene"); ax.set_ylabel("% of genes")
            continue
        kind, label, cap, bins = (("exons_per_transcript", "exons per transcript", 30, list(range(1, 32)))
                                  if panel == "exons" else
                                  ("cds_aa", "CDS length (aa)", 2000, list(range(0, 2050, 50))))
        for st_ in sets:
            vals = [min(v, cap) for v in series.get((st_, kind), [])]
            if vals:
                ax.hist(vals, bins=bins, histtype="step", lw=1.6, color=colors[st_], label=st_,
                        weights=[100 / len(vals)] * len(vals))
        ax.set_xlabel(f"{label}; values > {cap} shown at {cap}")
        ax.set_ylabel("% of transcripts")
    axes[0].legend(frameon=False, fontsize=8)
    f.tight_layout()
    iso_txt = (L("左：每个基因的转录本数（1、2、3、4、≥5 个）占基因的百分比；",
                 "Left: transcripts per gene (1, 2, 3, 4, ≥5) as % of genes; ") if iso_sets else
               L("（每个基因只有一个转录本，不画转录本数分布。）",
                 "(One transcript per gene: no transcripts-per-gene panel.) "))
    return fig_to_img(f, figdir, name, (
        L(f"{what}基因结构分布。", f"Gene structure distributions of the {what}."),
        iso_txt + L("中：每个转录本的外显子数；右：CDS 长度（氨基酸）。后两幅为各注释集内转录本的百分比，"
                    "便于比较规模不同的注释集；超过上限的值计入最后一格。",
                    "middle: exons per transcript; right: CDS length (amino acids). The last two show % of the "
                    "transcripts of each set, so that sets of different size can be compared; values above the limit "
                    "are counted in the last bin.")))


def ref_compare_help():
    if LANG == "zh":
        return """
<details open><summary><b>指标说明</b></summary>
<p>比较只看 CDS（编码区），不看 UTR：很多参考注释（尤其是 AUGUSTUS 等预测得到的）没有 UTR，比较 UTR 会把注释方式的差别误当成错误。
参考注释只取本次分析的序列，序列名按 <code>seqid_map.tsv</code> 换成流程内部名后再比较；所有比较都要求在同一条链上。</p>
<table><thead><tr><th>列</th><th>含义</th><th>怎么算</th></tr></thead><tbody>
<tr><td>ref_genes / pred_genes</td><td>参考 / 该注释集中的编码基因数</td><td>只计有 CDS 的基因</td></tr>
<tr><td>gene_recall_%</td><td>参考基因被找到的比例（宽松）</td><td>参考基因中，任一 CDS 片段与该注释集任一 CDS 有 ≥ 1 bp 同链重叠的比例</td></tr>
<tr><td>gene_precision_%</td><td>预测基因有参考对应的比例（宽松）</td><td>反过来：预测基因中，与参考 CDS 有同链重叠的比例。低值表示有较多参考里没有的基因（新基因、假阳性或参考漏注）</td></tr>
<tr><td>intron_recall_%</td><td>剪接位点的准确度</td><td>参考的 CDS 内含子（同一转录本相邻 CDS 片段之间的区间，按 序列+链+起点+终点 去重）中，被完全相同地预测出来的比例；差 1 bp 即不算</td></tr>
<tr><td>intron_precision_%</td><td>预测内含子的可靠程度</td><td>预测的 CDS 内含子中，与参考完全相同的比例</td></tr>
<tr><td>tx_recall_%</td><td>完整编码结构的准确度（最严格）</td><td>参考转录本的完整 CDS 结构（所有 CDS 片段的坐标，含起始和终止位置；按去重后的结构计）被完全相同地预测出来的比例；任何一个外显子边界不同都不算。单外显子基因也计入</td></tr>
<tr><td>tx_precision_%</td><td>预测的完整结构与参考一致的比例</td><td>预测转录本的完整 CDS 结构中，与参考完全相同的比例</td></tr>
</tbody></table>
<p><b>怎么读：</b>基因水平指标只要求有重叠，反映"基因找没找到"；内含子水平反映剪接位点是否准确；转录本水平要求整个编码结构完全一致，数值通常最低。
某个注释集召回高而精确低，说明它预测了较多参考里没有的基因或结构；反之说明它较保守。
主结果每个基因只保留一个转录本，参考注释若含多个异构体，其余异构体会拉低 tx_recall_%；相反，含全部异构体的注释集（如 refine）会拉低 tx_precision_%。
<b>参考注释本身也是预测结果时，这些数值反映的是两者的一致程度而不是准确率</b>，不一致处可能是任何一方的错误。</p>
</details>"""
    return """
<details open><summary><b>Metrics</b></summary>
<p>Only CDS (coding regions) are compared, not UTRs: many reference annotations (especially predictions such as AUGUSTUS) have no UTRs, and comparing UTRs would count differences in annotation style as errors.
Only the sequences of this analysis are taken from the reference, with sequence names converted to the internal names via <code>seqid_map.tsv</code>; every comparison requires the same strand.</p>
<table><thead><tr><th>Column</th><th>Meaning</th><th>How it is computed</th></tr></thead><tbody>
<tr><td>ref_genes / pred_genes</td><td>coding genes in the reference / in the annotation set</td><td>only genes with a CDS</td></tr>
<tr><td>gene_recall_%</td><td>share of reference genes found (lenient)</td><td>reference genes with any CDS segment overlapping any CDS of the set by ≥ 1 bp on the same strand</td></tr>
<tr><td>gene_precision_%</td><td>share of predicted genes with a reference counterpart (lenient)</td><td>the reverse: predicted genes overlapping a reference CDS on the same strand. Low values mean many genes absent from the reference (new genes, false positives or genes missing from the reference)</td></tr>
<tr><td>intron_recall_%</td><td>accuracy of splice sites</td><td>reference CDS introns (intervals between adjacent CDS segments of a transcript, unique by sequence + strand + start + end) predicted exactly; 1 bp off does not count</td></tr>
<tr><td>intron_precision_%</td><td>reliability of predicted introns</td><td>predicted CDS introns identical to a reference intron</td></tr>
<tr><td>tx_recall_%</td><td>accuracy of complete coding structures (strictest)</td><td>reference transcripts whose complete CDS structure (coordinates of all CDS segments, including start and stop; unique structures) is predicted exactly; any different exon border does not count. Single-exon genes are included</td></tr>
<tr><td>tx_precision_%</td><td>share of predicted complete structures identical to the reference</td><td>predicted transcripts whose complete CDS structure is identical to a reference one</td></tr>
</tbody></table>
<p><b>Reading:</b> gene-level metrics only require an overlap ("was the gene found"); intron level shows whether splice sites are exact; transcript level requires the whole coding structure to be identical and is usually lowest.
High recall with low precision means a set predicts many genes or structures absent from the reference; the reverse means it is conservative.
The main result keeps one transcript per gene: other isoforms of a reference with several isoforms lower tx_recall_%; conversely, a set with all isoforms (such as refine) lowers tx_precision_%.
<b>When the reference is itself a prediction, these values measure agreement, not accuracy</b>; a disagreement may be an error of either side.</p>
</details>"""


def tune_help():
    if LANG == "zh":
        return """
<details open><summary><b>原理</b></summary>
<p>EVM 按权重把各来源的证据（从头预测、同源投射、RNA 注释、蛋白比对、转录本组装）整合成一个共识模型，权重决定了冲突时更相信哪个来源。
配置中的默认权重是经验值，不一定适合当前物种和数据，所以流程在一条序列上用不同权重各跑一次 EVM，用原始证据给结果打分，选出最好的一组用于全基因组。</p>
<ol>
<li><b>调优区域：</b><code>evm.tune.region</code>，默认（auto）取最长的一条序列；EVM 的全部输入都截取到这条序列。</li>
<li><b>候选权重：</b>c00 为配置中的权重（基线）；其余为 <code>evm.tune.grid</code> 中各来源候选值的全部组合（未列在网格中的来源保持基线值，重复的组合去掉）。</li>
<li><b>打分</b>（每组权重的 EVM 结果）：
<ul>
<li><b>P（精确度）</b>：多外显子转录本中，CDS 内含子链被任一独立证据<b>完全</b>确认的比例。证据包括：长读（读段在模型范围内的内含子与模型完全一致，按全长 reads 加权，达到 <code>refine.min_lr_support</code>）；
没有长读时用短读（与修正步骤相同的标准：1 个内含子看剪接 reads，≥ 2 个内含子还要求 StringTie 组装出相同的链）；或一条与之完全一致的 miniprot 蛋白比对。
只要求"任一"证据，是为了不惩罚只有 RNA 证据、没有同源蛋白的物种特有基因。</li>
<li><b>R（召回）</b>：证据位点被 EVM 模型覆盖的比例，取下列可用项的平均：含剪接的长读中与某个 EVM CDS 重叠的比例；identity ≥ 0.8 的 miniprot 比对中与同链 EVM CDS 重叠的比例。</li>
<li><b>F = 2PR / (P + R)</b>，即 P 和 R 的调和平均。权重偏向某一来源时，可能多出很多基因（R 升高、P 下降）或丢掉基因（P 升高、R 下降），F 在两者之间取平衡。</li>
</ul></li>
<li><b>选择：</b>F 最高的一组；F 相同时取 P 较高者，再相同取编号较小者。<code>evm.tune.use_best: true</code> 时，全基因组 EVM 使用这组权重（见下方表格）；否则只报告结果，仍用配置中的权重。</li>
</ol>
<p><b>注意：</b>miniprot 比对和转录本组装本身也是 EVM 的输入，打分会略偏向给这两类证据较高权重的组合；只用一条序列打分，结果对全基因组有代表性但不完全等同。
各组 F 相差很小（例如 &lt; 0.005）时，可以认为这些权重效果相当，基线权重也是合理的选择。</p>
</details>"""
    return """
<details open><summary><b>Principle</b></summary>
<p>EVM combines the evidence of all sources (ab initio prediction, homology projection, RNA-based annotation, protein alignments, transcript assemblies) into consensus models according to weights; the weights decide which source is trusted more in a conflict.
The default weights in the configuration are rules of thumb and may not suit the species and data, so the workflow runs EVM once per weight set on one sequence, scores each result against the raw evidence, and uses the best set for the whole genome.</p>
<ol>
<li><b>Tuning region:</b> <code>evm.tune.region</code>, by default (auto) the longest sequence; all EVM inputs are cut to this sequence.</li>
<li><b>Candidate weights:</b> c00 = the weights of the configuration (baseline); the others are all combinations of the candidate values in <code>evm.tune.grid</code> (sources not in the grid keep their baseline value; duplicate combinations are removed).</li>
<li><b>Scores</b> (of the EVM result of each weight set):
<ul>
<li><b>P (precision)</b>: share of multi-exon transcripts whose CDS intron chain is <b>exactly</b> confirmed by any independent evidence: long reads (introns of the read within the model identical to the model, weighted by full-length reads, reaching <code>refine.min_lr_support</code>);
without long reads, short reads (same criteria as the refinement: junction reads for 1 intron, plus an identical StringTie chain for ≥ 2 introns); or an identical miniprot protein alignment.
"Any" evidence is required so that species-specific genes with RNA evidence but no homologous protein are not penalised.</li>
<li><b>R (recall)</b>: share of evidence sites covered by EVM models, the mean of the available items: spliced long reads overlapping an EVM CDS; miniprot alignments with identity ≥ 0.8 overlapping an EVM CDS on the same strand.</li>
<li><b>F = 2PR / (P + R)</b>, the harmonic mean of P and R. Weights favouring one source may add many genes (R up, P down) or lose genes (P up, R down); F balances the two.</li>
</ul></li>
<li><b>Choice:</b> the set with the highest F; ties go to the higher P, then to the lower number. With <code>evm.tune.use_best: true</code> the whole-genome EVM uses this set (table below); otherwise the result is only reported and the configured weights are used.</li>
</ol>
<p><b>Note:</b> miniprot alignments and transcript assemblies are EVM inputs themselves, so the scores slightly favour sets that weight these two higher; scoring on one sequence is representative of, but not identical to, the whole genome.
When the F values differ very little (e.g. &lt; 0.005), the weight sets perform equally and the baseline is a reasonable choice as well.</p>
</details>"""


REFINE_DEFAULTS = {"min_cds_aa": "50", "max_te_overlap": "0.5", "min_lr_support": "2", "single_exon_overlap": "0.5",
                   "sr_junction_fallback": "true", "min_sr_junction_reads": "3", "sr_min_junction_ratio": "0.1",
                   "sr_stringtie_min_cov": "3", "import_utr_isoforms": "true", "remove_unsupported_evm": "false"}


def config_section(config_path, section):
    """key -> value of a top-level YAML section, read without a YAML parser (report env has none)."""
    out, inside = {}, False
    if not config_path or not os.path.exists(config_path):
        return out
    for line in open(config_path):
        if re.match(r"^\S", line):
            inside = line.split(":")[0].strip() == section
            continue
        m = re.match(r"^  (\w+):\s*([^#\s][^#]*?)\s*(#.*)?$", line) if inside else None
        if m:
            out[m.group(1)] = m.group(2).strip().strip('"')
    return out


def refine_rules_html(config_path):
    """Collapsible full description of the refinement rules, with the values used in this run."""
    v = dict(REFINE_DEFAULTS)
    v.update({k: x for k, x in config_section(config_path, "refine").items() if k in REFINE_DEFAULTS})
    row = lambda *c: "<tr>" + "".join(f"<td>{x}</td>" for x in c) + "</tr>"
    unsup = (L("删除", "removed") if v["remove_unsupported_evm"] == "true"
             else L("保留并标记 refine_support=none", "kept and marked refine_support=none"))
    rules = "".join([
        row("R1", L("1 个 EVM 基因 : 1 个 EviAnn 基因，EviAnn 有一条转录本的 CDS 内含子链与 EVM 完全相同",
                    "1 EVM gene : 1 EviAnn gene, one EviAnn transcript has a CDS intron chain identical to EVM"),
            L("保留 EVM 的 CDS；外显子和 UTR 取自这条 EviAnn 转录本；EviAnn 的其它完整转录本（完整 ORF ≥ min_cds_aa）作为异构体加入"
              f"（import_utr_isoforms = {v['import_utr_isoforms']}）。若 EVM 的 CDS 不落在 EviAnn 外显子内，只保留 EVM。不依赖 RNA 支持",
              "EVM's CDS is kept; exons and UTRs come from that EviAnn transcript; the other complete EviAnn transcripts "
              f"(complete ORF ≥ min_cds_aa) are added as isoforms (import_utr_isoforms = {v['import_utr_isoforms']}). "
              "If EVM's CDS does not fall within the EviAnn exons, only EVM is kept. Does not depend on RNA support")),
        row("R2", L("只有 EviAnn 基因，没有重叠的 EVM 基因", "EviAnn gene only, no overlapping EVM gene"),
            L(f"同时满足才加入：① 有完整 ORF ≥ {v['min_cds_aa']} aa；② 有 RNA 证据（与 StringTie 转录本重叠，或"
              f"按下表“被支持”）；③ CDS 与 TE 重叠 ≤ {v['max_te_overlap']}。否则跳过（skip）",
              f"Added only when all hold: (1) a complete ORF ≥ {v['min_cds_aa']} aa; (2) RNA evidence (overlap with a "
              f"StringTie transcript, or \"supported\" as in the table below); (3) CDS overlap with TEs ≤ "
              f"{v['max_te_overlap']}. Otherwise skipped (skip)")),
        row("R3", L("1 : 1，但 CDS 内含子链不同", "1 : 1, but different CDS intron chains"),
            L("EviAnn 有被支持的完整转录本、EVM 不被支持 → 用 EviAnn 替换（replace_with_eviann）；两者都被支持 → 保留 EVM，"
              "并把被支持的 EviAnn 转录本加为异构体（keep_evm+eviann_isoforms）；其它 → 保留 EVM",
              "EviAnn has a supported complete transcript and EVM is not supported → replaced by EviAnn "
              "(replace_with_eviann); both supported → EVM kept and the supported EviAnn transcripts added as isoforms "
              "(keep_evm+eviann_isoforms); otherwise → EVM kept")),
        row(L("R4 拆分", "R4 split"), L("1 个 EVM 基因 : 多个 EviAnn 基因", "1 EVM gene : several EviAnn genes"),
            L("每个 EviAnn 基因都有被支持的转录本、且没有 reads 把它们连起来 → 用这些 EviAnn 基因替换；有 reads 连接 → 保留 EVM"
              "（说明确实是一个基因）；其它 → 保留 EVM，标记 refine_support=needs_review",
              "every EviAnn gene has a supported transcript and no reads connect them → replaced by these EviAnn genes; "
              "connecting reads → EVM kept (it is one gene); otherwise → EVM kept and marked refine_support=needs_review")),
        row(L("R4 融合", "R4 fusion"), L("多个 EVM 基因 : 1 个 EviAnn 基因", "several EVM genes : 1 EviAnn gene"),
            L("EviAnn 有被支持、且 CDS 覆盖所有这些 EVM 基因的转录本 → 用 EviAnn 替换；否则保留 EVM，标记 needs_review",
              "EviAnn has a supported transcript whose CDS covers all these EVM genes → replaced by EviAnn; otherwise EVM "
              "kept and marked needs_review")),
        row(L("R4 复杂", "R4 complex"), L("多 : 多", "many : many"),
            L("保留 EVM，标记 needs_review", "EVM kept, marked needs_review")),
        row("R5", L("只有 EVM 基因，没有重叠的 EviAnn 编码基因", "EVM gene only, no overlapping EviAnn coding gene"),
            L("有以下任一证据则保留并记录在 refine_support：蛋白（miniprot / GeMoMa 的 CDS 同链重叠）、RNA（StringTie 重叠或“被支持”）、"
              "BRAKER4 gene_support 的 RNA 支持；都没有 → ",
              "kept with any of these, recorded in refine_support: protein (CDS overlap with miniprot / GeMoMa on the "
              "same strand), RNA (StringTie overlap or \"supported\"), RNA support in BRAKER4 gene_support; none → ")
            + unsup + L(f"（remove_unsupported_evm = {v['remove_unsupported_evm']}）",
                        f" (remove_unsupported_evm = {v['remove_unsupported_evm']})")),
        row("R6", L("EviAnn 的 lncRNA", "EviAnn lncRNAs"),
            L("在更早的步骤已分出，单独输出为 lncRNA 文件，不参与修正",
              "separated in an earlier step, written to their own lncRNA file, not part of the refinement")),
    ])
    support = "".join([
        row(L("有长读", "long reads"), L("多外显子模型", "multi-exon model"),
            L("长读在模型范围内的内含子与模型的内含子完全一致；Iso-Seq HQ/LQ 转录本按 full_length_coverage（全长 reads 数）加权，"
              f"合计 ≥ min_lr_support = {v['min_lr_support']}",
              "introns of the long read within the model identical to the model's introns; Iso-Seq HQ/LQ transcripts "
              f"weighted by full_length_coverage (full-length reads), total ≥ min_lr_support = {v['min_lr_support']}")),
        row(L("有长读", "long reads"), L("单外显子模型", "single-exon model"),
            L(f"无剪接的长读覆盖模型长度 ≥ single_exon_overlap = {v['single_exon_overlap']}（加权计数同上）",
              f"unspliced long reads cover ≥ single_exon_overlap = {v['single_exon_overlap']} of the model "
              "(weighted as above)")),
        row(L("只有短读（B）", "short reads only (B)"), L("1 个内含子", "1 intron"),
            L(f"该内含子的剪接 reads ≥ min_sr_junction_reads = {v['min_sr_junction_reads']}，且 ≥ sr_min_junction_ratio = "
              f"{v['sr_min_junction_ratio']} × 与它重叠的最强竞争内含子的 reads 数（排除低频的错误剪接）",
              f"junction reads of the intron ≥ min_sr_junction_reads = {v['min_sr_junction_reads']}, and ≥ "
              f"sr_min_junction_ratio = {v['sr_min_junction_ratio']} x the reads of the strongest overlapping "
              "competing intron (excludes rare mis-splicing)")),
        row(L("只有短读（A 且 B）", "short reads only (A and B)"), L("≥ 2 个内含子", "≥ 2 introns"),
            L("每个内含子都满足 B；并且有一条 StringTie 转录本（短读组装，cov ≥ sr_stringtie_min_cov = "
              f"{v['sr_stringtie_min_cov']}，同链或无链）在模型范围内的内含子链与模型完全一致，证明这些内含子属于同一条转录本",
              "every intron meets B; and a StringTie transcript (short-read assembly, cov ≥ sr_stringtie_min_cov = "
              f"{v['sr_stringtie_min_cov']}, same strand or unstranded) has an intron chain within the model identical "
              "to the model, showing that the introns belong to one transcript")),
        row(L("只有短读（C）", "short reads only (C)"), L("单外显子模型", "single-exon model"),
            L(f"有单外显子 StringTie 转录本（cov 达标）覆盖模型长度 ≥ {v['single_exon_overlap']}，且模型范围内没有达到 B 绝对门槛的剪接位点",
              f"a single-exon StringTie transcript (sufficient cov) covers ≥ {v['single_exon_overlap']} of the model, "
              "and no splice site within the model reaches the absolute threshold of B")),
        row(L("无 RNA 数据", "no RNA data"), "—", L("一律视为不被支持", "always counted as not supported")),
    ])
    th = lambda *c: "<thead><tr>" + "".join(f"<th>{x}</th>" for x in c) + "</tr></thead>"
    return ("<details><summary><b>" + L("修正规则详解（点击展开）", "Refinement rules in detail (click to expand)") +
            "</b></summary><p><b>" + L("位点与比较方式：", "Loci and comparison: ") + "</b>" +
            L("EVM 基因和 EviAnn 编码基因若在同一条链上 CDS 有重叠，归入同一个位点（可连锁合并），"
              "规则按位点内两边的基因数选择。模型之间比较 CDS 的内含子链（EVM 模型没有 UTR，比较完整外显子很难一致）。"
              f"EviAnn 转录本须有完整 ORF（起始和终止密码子、无内部终止）且 ≥ min_cds_aa = {v['min_cds_aa']} aa 才会被采用。",
              "EVM genes and EviAnn coding genes whose CDS overlap on the same strand form one locus (chained), and the "
              "rule is chosen by the number of genes on each side. Models are compared by their CDS intron chains (EVM "
              "models have no UTRs, so complete exons rarely agree). An EviAnn transcript is used only with a complete "
              f"ORF (start and stop codon, no internal stop) of ≥ min_cds_aa = {v['min_cds_aa']} aa.") + "</p>"
            "<h4>" + L("规则", "Rules") + "</h4><table>" + th(L("规则", "Rule"), L("情形", "Case"), L("处理", "Action")) +
            "<tbody>" + rules + "</tbody></table>"
            "<h4>" + L("“被支持”的判定（RNA 支持）", "When a model is \"supported\" (RNA support)") + "</h4><table>" +
            th(L("数据", "Data"), L("模型", "Model"), L("条件", "Condition")) + "<tbody>" + support + "</tbody></table>"
            "<p><b>" + L("跨基因连接（R4）：", "Connections between genes (R4): ") + "</b>" +
            L("有长读时为同时含两个相邻基因外显子的长读；只有短读时为一端在基因 A 外显子、另一端在基因 B 外显子的剪接 reads。",
              "with long reads, a long read containing exons of both neighbouring genes; with short reads only, junction "
              "reads with one end in an exon of gene A and the other in an exon of gene B.") + "</p>"
            "<p><b>" + L("输出记录：", "Records: ") + "</b>" +
            L("基因属性 refine_rule（适用的规则）、eviann_support（是否与 EviAnn 一致或来自 EviAnn）、refine_support"
              "（R5 的证据或 needs_review）；05_EVM/refine/refine_decisions.tsv 每个位点一行，support_basis 列记录支持依据（long / A+B / B / C）。"
              "表中数值为本次运行的配置（refine 段）。",
              "gene attributes refine_rule (the rule applied), eviann_support (agrees with or comes from EviAnn), "
              "refine_support (R5 evidence or needs_review); 05_EVM/refine/refine_decisions.tsv has one line per locus, "
              "the support_basis column gives the basis of the support (long / A+B / B / C). The values in the tables "
              "are this run's configuration (refine section).") + "</p></details>")


def long_read_weights_html(path):
    """Full-length read weights of the long-read samples (refine / tuning support), with a warning when Iso-Seq
    data carry no weights (HQ transcripts whose headers lost full_length_coverage, and no count table)."""
    h, r = read_tsv(path)
    if not h:
        return ""
    i_p, i_w = h.index("platform"), h.index("weighted_reads")
    warn = [x[0] for x in r if x[i_p] == "pacbio_isoseq" and x[i_w] == "0"]
    note = ""
    if warn:
        names = html.escape(L("、", ", ").join(warn))
        note = ("<div class='note'><b>" + L("注意：", "Note: ") + "</b>" +
                L(f"Iso-Seq 样本 {names} 没有任何一条记录带全长 reads 数"
                  "（表头没有 full_length_coverage，样本表也没有 counts 列），每条记录按 1 计数。若这些记录是原始 reads（CCS / FLNC），"
                  "这是正确的；若是聚类后的 HQ 转录本，长读支持会被明显低估（修正 R3/R4 偏保守、调优 P 偏低）。"
                  "请在 rna_long.tsv 的 counts 列填写 isoseq cluster 的 cluster_report.csv 或 isoseq collapse 的 flnc_count.txt，"
                  "或将 refine.min_lr_support 设为 1。",
                  f"no record of the Iso-Seq sample(s) {names} carries a full-length read count (no "
                  "full_length_coverage in the headers and no counts column in the sample table), so every record "
                  "counts 1. This is correct for raw reads (CCS / FLNC); for clustered HQ transcripts the long-read "
                  "support is clearly underestimated (refinement R3/R4 too conservative, tuning P too low). Fill the "
                  "counts column of rna_long.tsv with isoseq cluster's cluster_report.csv or isoseq collapse's "
                  "flnc_count.txt, or set refine.min_lr_support to 1.") + "</div>")
    return ("<h3>" + L("长读的全长 reads 权重", "Full-length read weights of long reads") + "</h3><p class='muted'>" +
            L("修正步骤（R3/R4）和权重调优按\"全长 reads\"计算长读支持："
              "每条比对上的记录的权重取自样本表 counts 列的计数文件，或 FASTA/FASTQ 表头的 full_length_coverage，都没有时按 1 计。",
              "The refinement (R3/R4) and the weight tuning count long-read support in full-length reads: the weight "
              "of each aligned record comes from the count file in the counts column of the sample table, or from "
              "full_length_coverage in the FASTA/FASTQ header, and is 1 otherwise.") + "</p>" + note + table(h, r))


def tune_html(path, best, figdir):
    """EVM weight-tuning results: score table, F bar chart, the weights finally used."""
    if not os.path.exists(path):
        return "<p class='muted'>" + L("未进行权重调优（evm.tune.enabled: false），使用配置中的权重。",
                                       "No weight tuning (evm.tune.enabled: false); the configured weights are used.") + "</p>"
    lines = [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]
    note = [l[0].lstrip("# ") for l in lines if l[0].startswith("#")]
    body = [l for l in lines if not l[0].startswith("#")]
    if len(body) < 2:
        return ""
    h, r = body[0], body[1:]
    fig = ""
    if plt is not None:
        ix, ip, ir = h.index("F"), h.index("P"), h.index("R")
        rows_ = r[:20]
        f, ax = plt.subplots(figsize=(8.5, 0.42 * len(rows_) + 1.3))
        y = list(range(len(rows_)))
        for col, mk, colr, lab in ((ip, "o", "#e0a43a", "P"), (ir, "s", "#5a9e6f", "R"), (ix, "D", "#2f6fb0", "F")):
            ax.scatter([float(x[col]) for x in rows_], y, marker=mk, color=colr, label=lab, zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{x[1]}  {x[h.index('weights')]}" for x in rows_], fontsize=8)
        ax.invert_yaxis(); ax.grid(axis="x", color="#ddd", lw=0.6)
        ax.set_xlabel("score")
        legend_below(ax, ncol=3)
        fig = fig_to_img(f, figdir, "evm_tuning", (
            L("EVM 权重调优。", "EVM weight tuning."),
            L("每行一组权重（c00 为配置中的基线），按 F 从高到低排列；P = 多外显子模型的内含子链被任一证据完全确认的比例，"
              "R = 证据位点被 EVM 模型覆盖的比例，F = 2PR/(P+R)。横轴只显示数值所在区间，便于看出差别。",
              "One row per weight set (c00 = the configured baseline), sorted by F; P = share of multi-exon models "
              "whose intron chain is exactly confirmed by any evidence, R = share of evidence sites covered by EVM "
              "models, F = 2PR/(P+R). The x axis only spans the range of the values, so that differences are visible.")))
    used = [l.rstrip("\n").split("\t") for l in open(best) if l.strip()] if os.path.exists(best) else []
    base = next((x for x in r if x[1] == "c00"), None)
    cmp_ = (L(f"最优组合 <b>{html.escape(r[0][1])}</b>（F = {r[0][h.index('F')]}）；配置中的基线权重 c00 排第 "
              f"{base[0]} 位（F = {base[h.index('F')]}）。",
              f"Best set <b>{html.escape(r[0][1])}</b> (F = {r[0][h.index('F')]}); the configured baseline c00 ranks "
              f"{base[0]} (F = {base[h.index('F')]}).") if base else "")
    return ("<h3>" + L("权重调优", "Weight tuning") + "</h3>" + tune_help() + "<p class='muted'>" +
            html.escape("; ".join(note)) + "</p><div class='note'>" + cmp_ + "</div>" + fig + table(h, r, limit=40) +
            "<p>" + L("全基因组 EVM 实际使用的权重：", "Weights used by the whole-genome EVM:") + "</p>" +
            table([L("类别", "Class"), L("来源", "Source"), L("权重", "Weight")], used))


def naming_html(id_map, prefix):
    """Genes per chromosome code, from the old -> new ID table."""
    if not os.path.exists(id_map):
        return ""
    per, first = {}, {}
    with open(id_map) as fh:
        next(fh)
        for line in fh:
            kind, old, new, seqid = line.rstrip("\n").split("\t")[:4]
            if kind != "gene" or not new.startswith(prefix):
                continue
            m = re.match(r"(.*?)([GL])\d+$", new[len(prefix):])
            if not m:
                continue
            code, tag = m.group(1) or "-", m.group(2)
            per.setdefault(code, {"G": 0, "L": 0, "seqs": set()})
            per[code][tag] += 1
            per[code]["seqs"].add(seqid)
            first.setdefault((code, tag), new)
    if not per:
        return ""
    rows = []
    for code in sorted(per, key=lambda c: (c == "00", c)):
        v = per[code]
        seqs = sorted(v["seqs"])
        more = L(f" 等 {len(seqs)} 条", f" ... ({len(seqs)} in total)") if len(seqs) > 3 else ""
        rows.append([code, ", ".join(seqs[:3]) + more, v["G"], v["L"], first.get((code, "G"), "-")])
    return ("<h3>" + L("基因命名与各染色体基因数", "Gene names and genes per chromosome") + "</h3><p class='muted'>" +
            L("编码基因 &lt;前缀&gt;&lt;染色体号&gt;G&lt;编号&gt;，"
              "lncRNA 用 L；未挂载到染色体的序列代码为 00；新旧 ID 对照见 05_EVM/final/id_map.tsv。",
              "Coding genes &lt;prefix&gt;&lt;chromosome&gt;G&lt;number&gt;, lncRNAs use L; sequences not assigned "
              "to a chromosome have code 00; old and new IDs are listed in 05_EVM/final/id_map.tsv.") + "</p>" +
            table([L("染色体代码", "Chromosome code"), L("序列", "Sequences"), L("编码基因", "Coding genes"), "lncRNA",
                   L("第一个基因", "First gene")], rows))


def source_support_html(support_tsv, combos_tsv, figdir, name):
    sh_, sr_ = read_tsv(support_tsv)
    ch_, cr_ = read_tsv(combos_tsv)
    if not sh_:
        return ""
    fig = ""
    if cr_:
        top = cr_[:12]
        what = (L("EVM 结果", "EVM result") if name == "evm" else L("最终结果", "final result"))
        fig = bar_figure([x[0] for x in top], [int(x[1]) for x in top], "transcripts", figdir,
                         f"{name}_source_combinations", (
                             L(f"{what}的来源组合。", f"Source combinations of the {what}."),
                             L("每个转录本按哪些输入来源（BRAKER4 / GeMoMa / EviAnn）有 CDS 结构完全相同的模型分组；"
                               "no exact match = 与任何来源都不完全一致（EVM 拼合或修正产生的结构）。",
                               "Transcripts grouped by the input sources (BRAKER4 / GeMoMa / EviAnn) that have a model "
                               "with an identical CDS structure; no exact match = identical to no source (structures "
                               "assembled by EVM or produced by the refinement).")), values_text="{:.0f}")
    return ("<h3>" + L("与各输入来源的一致性", "Agreement with the input sources") + "</h3><p class='muted'>" +
            L("exact = 该来源有 CDS 结构完全相同的转录本；overlap = 与该来源的 CDS 在同一条链上有重叠。组合 = 有完全一致模型的来源组合。",
              "exact = the source has a transcript with an identical CDS structure; overlap = the CDS overlaps one of "
              "the source on the same strand. Combination = the sources with an identical model.") + "</p>" +
            table(sh_, sr_) + fig + table(ch_, cr_, limit=30))


def section_what(name):
    return {"eviann": L("EviAnn 注释（编码基因与 lncRNA）", "EviAnn annotation (coding genes and lncRNAs)"),
            "braker4": L("BRAKER4 注释", "BRAKER4 annotation"), "evm": L("EVM 结果", "EVM result"),
            "gemoma": L("同源注释（各物种单独结果与合并结果）", "homology annotation (per species and merged)"),
            "final": L("最终结果（主结果）", "final result (main result)")}.get(name, "")


def annotation_section(O, name, stats_tsv, dist_tsv, qc_dir, figdir, intro, extra=""):
    """Statistics table + distributions + protein-set quality + AGAT for one annotation source."""
    h, r = read_tsv(stats_tsv)
    _, dist = read_tsv(dist_tsv)
    parts = [intro, "<h3>" + L("基因结构统计", "Gene structure statistics") + "</h3>",
             "<p class='muted'>" + L("UTR 优先取 GFF 中的 UTR 特征，否则按外显子超出 CDS 的部分判断（区分正负链）。",
                                     "UTRs are taken from UTR features of the GFF when present, otherwise from the "
                                     "parts of exons outside the CDS (strand-aware).") + "</p>",
             table(h, r), distribution_figure(dist, figdir, f"{name}_distributions", section_what(name))]
    q = protein_qc_rows(qc_dir)
    if q:
        parts.append("<h3>" + L("蛋白集质量（最长异构体）", "Protein set quality (longest isoform)") + "</h3>" +
                     table([L("评估", "Assessment"), L("结果", "Result")], q))
    parts.append(extra)
    agat_f = os.path.join(qc_dir, "agat", "stats.txt")
    if os.path.exists(agat_f):
        parts.append("<details><summary>" + L("AGAT 统计全文", "AGAT statistics (full text)") + "</summary><pre>" +
                     html.escape(open(agat_f).read()[:20000]) + "</pre></details>")
    return "".join(parts)


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _hms(sec):
    if sec is None:
        return "-"
    sec = int(round(sec))
    return f"{sec // 3600}:{sec % 3600 // 60:02d}:{sec % 60:02d}"


def _gb(mb):
    return "-" if mb is None else f"{mb / 1024:.1f}"


CONTAINER_RULES = {"braker4"}     # tools run in Apptainer containers: the benchmark misses their CPU / memory


def resources_html(outdir, requests_tsv, figdir):
    """Runtime and resources of every job (Snakemake benchmark files <step>/logs/benchmarks/*.tsv) against the
    requested threads / memory / runtime (07_Report/resource_requests.tsv), per rule with suggestions."""
    files = sorted(glob.glob(os.path.join(outdir, "*", "logs", "benchmarks", "*.tsv")))
    if not files:
        return ""
    rh, rr = read_tsv(requests_tsv)
    pats = []                                    # (regex of the benchmark path, rule, threads, mem_mb, runtime)
    if rh:
        for x in rr:
            rx = re.escape(x[4])
            rx = re.sub(r"\\\{[^}]*\\\}", "[^/]+", rx)  # escaped {wildcard} -> one path component part
            pats.append((re.compile(rx + "$"), x[0], _num(x[1]), _num(x[2]), _num(x[3])))
    jobs = []
    for f_ in files:
        h, r = read_tsv(f_)
        if not h or not r:
            continue
        v = dict(zip(h, r[-1]))                  # last repeat
        rel = os.path.relpath(f_, outdir)
        rule, thr, mem, rt = next(((p[1], p[2], p[3], p[4]) for p in pats if p[0].match(rel)), ("?", None, None, None))
        wall, cpu, rss = _num(v.get("s")), _num(v.get("cpu_time")), _num(v.get("max_rss"))
        eff = (100 * cpu / (wall * thr)) if (cpu is not None and wall and thr) else None
        jobs.append({"step": rel.split("/")[0], "job": os.path.basename(f_)[:-4], "rule": rule, "threads": thr,
                     "mem": mem, "rt": rt, "wall": wall, "cpu": cpu, "rss": rss, "eff": eff})
    if not jobs:
        return ""
    total_wall = sum(j["wall"] or 0 for j in jobs)
    # per rule: worst case over its jobs
    per = {}
    for j in jobs:
        per.setdefault(j["rule"], []).append(j)
    rows, sugg_n = [], 0
    for rule, js in sorted(per.items(), key=lambda t: -max(j["wall"] or 0 for j in t[1])):
        thr, mem, rt = js[0]["threads"], js[0]["mem"], js[0]["rt"]
        wmax = max((j["wall"] or 0) for j in js)
        rss = max((j["rss"] or 0) for j in js) or None
        cores = max(((j["cpu"] or 0) / j["wall"]) for j in js if j["wall"]) if any(j["wall"] for j in js) else None
        eff = (100 * cores / thr) if (cores is not None and thr) else None
        tips = []
        if rss and mem:
            if rss > 0.9 * mem:
                tips.append(L(f"内存接近或超过申请值，建议 ≥ {int(rss * 1.3 / 1024) + 1}G",
                              f"memory near or above the request: use ≥ {int(rss * 1.3 / 1024) + 1}G"))
            elif rss < 0.4 * mem and mem >= 4096:
                tips.append(L(f"内存用量远低于申请值，可降到约 {max(1, int(rss * 1.5 / 1024) + 1)}G",
                              f"memory far below the request: about {max(1, int(rss * 1.5 / 1024) + 1)}G is enough"))
        if eff is not None and thr and thr >= 4 and wmax >= 300 and eff < 40:
            tips.append(L(f"平均只用到约 {cores:.1f} 个核，可减少 threads",
                          f"uses about {cores:.1f} cores on average: fewer threads would do"))
        if rt and wmax > 0.8 * rt * 60:
            tips.append(L("运行时间接近申请值，建议增大 runtime", "run time close to the request: raise runtime"))
        if rule in CONTAINER_RULES:                 # CPU / memory of the container processes are not measured
            tips = [t for t in tips if "runtime" in t] + [
                L("容器内的进程未计入 CPU 和内存（实测 BRAKER4 峰值可达 100 GB 级），请以 SLURM 的 sacct 为准",
                  "CPU and memory of the processes in the containers are not counted (BRAKER4 peaks were ~100 GB in a "
                  "server test): use SLURM's sacct instead")]
        sugg_n += bool(tips)
        rows.append([rule, len(js), _hms(wmax), _hms(sum(j["wall"] or 0 for j in js)),
                     "-" if thr is None else int(thr), "-" if cores is None else f"{cores:.1f}",
                     "-" if eff is None else f"{eff:.0f}", _gb(rss), _gb(mem) if mem else L("默认", "default"),
                     "-" if (not rss or not mem) else f"{100 * rss / mem:.0f}",
                     _hms(rt * 60) if rt else "-", L("；", "; ").join(tips)])
    head = [L("规则", "Rule"), L("作业数", "Jobs"), L("最长用时", "Longest"), L("合计用时", "Total"),
            L("申请线程", "Threads"), L("实际核数", "Cores used"), L("CPU 效率 %", "CPU eff. %"),
            L("峰值内存 GB", "Peak mem GB"), L("申请内存 GB", "Mem req. GB"), L("内存占比 %", "Mem used %"),
            L("申请时间", "Runtime req."), L("建议", "Suggestion")]
    top = sorted(jobs, key=lambda j: -(j["wall"] or 0))[:20]
    fig = bar_figure([f"{j['step']}/{j['job']}" for j in top], [(j["wall"] or 0) / 3600 for j in top], "hours",
                     figdir, "runtime_top", (
                         L("用时最长的步骤。", "Longest-running steps."),
                         L("实际运行时间（墙钟时间，小时）最长的 20 个作业。", "The 20 jobs with the longest wall-clock "
                                                                          "time (hours).")), wrap=60,
                     values_text="{:.2f}")
    detail = [[j["step"], j["job"], j["rule"], _hms(j["wall"]), _hms(j["cpu"]),
               "-" if j["threads"] is None else int(j["threads"]), "-" if j["eff"] is None else f"{j['eff']:.0f}",
               _gb(j["rss"]), _gb(j["mem"]) if j["mem"] else L("默认", "default")]
              for j in sorted(jobs, key=lambda j: -(j["wall"] or 0))]
    return ("<p>" + L("每个作业的实际运行时间、CPU 时间和峰值内存由 Snakemake benchmark 记录（"
                      "&lt;步骤&gt;/logs/benchmarks/*.tsv），与申请的资源（config.yaml 的 resources，未设置的用 profile 默认值）对比，"
                      "供下次调整 resources 参考。",
                      "Wall-clock time, CPU time and peak memory of every job are recorded by Snakemake benchmarks "
                      "(&lt;step&gt;/logs/benchmarks/*.tsv) and compared with the requested resources (resources in "
                      "config.yaml; profile defaults where not set), as a guide for adjusting resources next time.") +
            "</p><div class='note'>" +
            L(f"共 {len(jobs)} 个作业，墙钟时间合计 {_hms(total_wall)}（各作业相加，并行运行时实际耗时更短）；"
              f"{sugg_n} 条规则有调整建议。",
              f"{len(jobs)} jobs, {_hms(total_wall)} wall-clock time in total (summed over jobs; less in practice when "
              f"they run in parallel); suggestions for {sugg_n} rules.") + "</div>"
            "<p class='muted'>" +
            L("实际核数 = CPU 时间 / 墙钟时间；CPU 效率 = 实际核数 / 申请线程。峰值内存为 max_rss（含子进程）。"
              "建议按本次运行的最大值给出（内存留 30–50% 余量）；数据规模不同（例如测试用单条染色体，正式用全基因组）时需按比例放大。"
              "单机运行时内存和时间只用于调度，不会被强制限制；SLURM 上超过申请值的作业会被终止。"
              "BRAKER4 以 executor: slurm 运行时，其子作业不计入。",
              "Cores used = CPU time / wall-clock time; CPU efficiency = cores used / requested threads. Peak memory is "
              "max_rss (including child processes). Suggestions are based on the maxima of this run (memory with "
              "30–50% headroom); scale them for other data sizes (e.g. a single chromosome in a test versus the whole "
              "genome). Running locally, memory and time are only used for scheduling and not enforced; on SLURM jobs "
              "exceeding the request are killed. With braker4.executor: slurm, BRAKER4's own jobs are not counted.") +
            "</p><h3>" + L("按规则汇总", "Per rule") + "</h3>" + table(head, rows) + fig +
            "<details><summary>" + L("全部作业明细", "All jobs") + "</summary>" +
            table([L("步骤", "Step"), L("作业", "Job"), L("规则", "Rule"), L("墙钟时间", "Wall time"),
                   L("CPU 时间", "CPU time"), L("申请线程", "Threads"), L("CPU 效率 %", "CPU eff. %"),
                   L("峰值内存 GB", "Peak mem GB"), L("申请内存 GB", "Mem req. GB")], detail, limit=1000) +
            "</details>")


EVIANN_PROTEINS = {"swissprot": ("UniProt-SwissProt（eviann.swissprot）", "UniProt-SwissProt (eviann.swissprot)"),
                   "homologs": ("同源物种蛋白（未设置 eviann.swissprot）",
                                "homologous species proteins (eviann.swissprot not set)"),
                   "none": ("未运行（无 RNA）", "not run (no RNA)")}


def build(a, out):
    """Build the report in the current LANG and write it to `out`."""
    O = lambda *p: os.path.join(a.outdir, *p)
    figdir = os.path.join(os.path.dirname(out), "figures")
    S = []

    # 00_Data: inputs
    rn = kv_file(O("00_Data/genome/rename_report.tsv"))
    hom = sorted(glob.glob(O("00_Data/qc/homologs/*.clean_report.tsv")))
    hom_html = ""
    if hom:
        items, per = [], {}
        for f_ in hom:
            name = os.path.basename(f_)[:-len(".clean_report.tsv")]
            per[name] = dict(kv_file(f_))
            for k in per[name]:
                if k not in items and k != "item":
                    items.append(k)
        hom_html = ("<h3>" + L("同源物种注释预处理", "Preprocessing of the homologous species' annotations") +
                    "</h3><p class='muted'>" +
                    L("只保留蛋白编码基因（及其有 CDS 的转录本），"
                      "去掉转座子基因、假基因、ncRNA 基因等；GeMoMa 和蛋白库都使用预处理后的注释。",
                      "Only protein-coding genes (and their transcripts with a CDS) are kept; transposon genes, "
                      "pseudogenes, ncRNA genes etc. are removed; GeMoMa and the protein database use the "
                      "preprocessed annotations.") + "</p>" +
                    table([L("项目", "Item")] + list(per), [[k] + [per[n].get(k, "0") for n in per] for k in items]))
    S.append(("data", L("00_Data 输入与序列名检查", "00_Data Inputs and sequence names"),
              table([L("项目", "Item"), L("值", "Value")], rn) + hom_html))

    # 00_Data: genome assembly quality (every assessment result appears here)
    AQ = lambda *p: O("00_Data/qc/assembly", *p)
    qh, qr = read_tsv(O("00_Data/qc/assembly_quality.tsv"))
    parts = ["<p>" + L("对用户提供的原始组装（原始序列名）进行评估。连续性：scaffold / contig 两级，contig 在连续 N 处拆分；"
                       "完整度：BUSCO 基因组模式、Merqury k-mer 完整度（可选）、LAI（EDTA 完成后计算）；准确度：Merqury QV（可选）；"
                       "端粒：tidk；污染：FCS-GX（可选，只报告不删除）。",
                       "The assembly as provided (original sequence names) is assessed. Contiguity: scaffold and "
                       "contig level, contigs split at runs of N; completeness: BUSCO genome mode, Merqury k-mer "
                       "completeness (optional), LAI (after EDTA); accuracy: Merqury QV (optional); telomeres: tidk; "
                       "contamination: FCS-GX (optional, reported only, nothing removed).") + "</p>",
             "<h3>" + L("汇总", "Summary") + "</h3>" + table(qh, qr)]
    skipped = [x[1] for x in qr if len(x) > 2 and x[2] == "not run"]
    if skipped:
        parts.append("<div class='note'>" + L("未运行（未配置所需数据或数据库）：", "Not run (required data or database "
                                                                         "not configured): ") +
                     html.escape(L("、", ", ").join(skipped)) + "</div>")

    sh_, sr_ = read_tsv(AQ("assembly_stats.tsv"))
    _, seqs = read_tsv(AQ("sequences.tsv"))
    fig = ""
    if plt is not None and len(seqs) >= 2:          # one sequence: the curve is a single point
        lens = sorted((int(x[1]) for x in seqs), reverse=True)
        total, acc, xs, ys = sum(lens), 0, [], []
        for i, l in enumerate(lens, 1):
            acc += l
            xs.append(i)
            ys.append(100 * acc / total)
        f, ax = plt.subplots(figsize=(7, 3.5))
        ax.plot(xs, ys, color="#2f6fb0")
        ax.axhline(50, color="#999", lw=0.8, ls="--")
        ax.axhline(90, color="#999", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel("sequences (longest first, log scale)")
        ax.set_ylabel("cumulative % of assembly")
        fig = fig_to_img(f, figdir, "assembly_cumulative_length", (
            L("序列累计长度。", "Cumulative sequence length."),
            L("序列按长度从长到短排列（横轴为对数刻度），纵轴为累计长度占组装总长的百分比；"
              "虚线 50% 和点线 90% 与曲线的交点分别对应 L50 和 L90。",
              "Sequences sorted from longest to shortest (log x axis); the y axis is the cumulative length as % of the "
              "assembly; the curve crosses the dashed 50% and dotted 90% lines at L50 and L90.")))
    sk_h, sk_r = read_tsv(O("00_Data/qc/seqkit/genome.stats.tsv"))
    parts.append("<h3>" + L("连续性", "Contiguity") + "</h3>" + table(sh_, sr_) + fig +
                 "<p class='muted'>" + L("seqkit stats -a（重命名和长度过滤后的内部基因组）：",
                                         "seqkit stats -a (internal genome after renaming and length filter):") +
                 "</p>" + table(sk_h, sk_r))
    if seqs:
        parts.append("<h3>" + L("各序列", "Sequences") + "</h3>" +
                     table(["seq_id", "length", "GC_percent", "N_bases", "gaps", "gap_length", "contigs"], seqs, limit=60))

    # BUSCO genome mode
    brows = []
    for f_ in sorted(glob.glob(O("00_Data/qc/busco_genome/*/short_summary.txt"))):
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", open(f_).read())
        if m:
            brows.append(["genome", os.path.basename(os.path.dirname(f_)), *m.groups()])
    if brows:
        bh = ["source", "lineage", "C", "S", "D", "F", "M", "n"]
        parts.append("<h3>" + L("BUSCO（基因组模式）", "BUSCO (genome mode)") + "</h3>" +
                     busco_figure(bh, brows, figdir, "busco_genome", (
                         L("基因组 BUSCO（genome 模式）。", "Genome BUSCO (genome mode)."),
                         L("每行一个 BUSCO 数据集；按单拷贝完整（S）、重复完整（D）、片段（F）、缺失（M）的百分比堆叠，"
                           "右侧为完整度 C = S + D。这是基因组组装完整性的基线，注释的 BUSCO 一般不会高于它。",
                           "One row per BUSCO dataset; bars stack complete single-copy (S), complete duplicated (D), "
                           "fragmented (F) and missing (M); completeness C = S + D on the right. This is the baseline "
                           "of assembly completeness; the BUSCO of an annotation is usually not higher."))) +
                     table(bh, brows))

    # Merqury
    mrows = []
    for name, cols in (("genome.qv", ["assembly", "kmers_asm_only", "kmers_total", "QV", "error_rate"]),
                       ("genome.completeness.stats", ["assembly", "set", "kmers_asm", "kmers_total", "completeness_%"])):
        f_ = AQ("merqury", name)
        if os.path.exists(f_):
            r = [l.split() for l in open(f_) if l.strip()]
            mrows.append((cols, r))
    if mrows:
        parts.append("<h3>Merqury</h3>" + "".join(table(c, r) for c, r in mrows))

    # LAI (after EDTA)
    lai_f = AQ("lai", "genome.LAI")
    if os.path.exists(lai_f):
        lai_help = ("<p class='muted'>" +
                    L("LAI（LTR Assembly Index，Ou et al. 2018, NAR 46:e126）用完整 LTR 反转录转座子占全部 LTR 序列的比例"
                      "衡量组装对重复区的连续性：raw_LAI = 完整 LTR 长度 / 全部 LTR 长度 × 100；LAI 再按 LTR 序列一致性校正，"
                      "排除近期大量扩增的家族造成的偏高。参考分级：draft &lt; 10，reference 10–20，gold ≥ 20。"
                      "完整 LTR 来自 EDTA，全部 LTR 来自 RepeatMasker（EDTA 库）。",
                      "LAI (LTR Assembly Index, Ou et al. 2018, NAR 46:e126) measures the contiguity of repeat regions "
                      "as the share of intact LTR retrotransposons in all LTR sequence: raw_LAI = intact LTR length / "
                      "total LTR length x 100; LAI is further corrected for LTR identity to exclude inflation by recently "
                      "amplified families. Guide values: draft &lt; 10, reference 10–20, gold ≥ 20. Intact LTRs come "
                      "from EDTA, all LTRs from RepeatMasker (EDTA library).") + "</p>")
        first = open(lai_f).readline()
        if first.startswith("# LAI not applicable"):
            reason = html.escape(first.split(":", 1)[1].strip().rstrip("."))
            parts.append("<h3>LAI</h3>" + lai_help + "<div class='note'>" +
                         L(f"LAI 不适用：{reason}。LTR 反转录转座子占基因组不足 5% 时，LAI 无法可靠估计（单条染色体测试或 LTR 很少的物种常见）。",
                           f"LAI not applicable: {reason}. With LTR retrotransposons below 5% of the genome LAI cannot "
                           "be estimated reliably (common for single-chromosome tests or species with few LTRs).") +
                         "</div>")
        else:
            lh, lr = read_tsv_ws(lai_f)
            whole = [x for x in lr if x and x[0] == "whole_genome"]
            parts.append("<h3>LAI</h3>" + lai_help +
                         ("<div class='note'>" + L("全基因组 LAI = ", "Whole-genome LAI = ") +
                          f"<b>{html.escape(whole[0][-1])}</b>" +
                          L(f"（raw_LAI = {html.escape(whole[0][-2])}）", f" (raw_LAI = {html.escape(whole[0][-2])})") +
                          "</div>" if whole else "") +
                         "<details><summary>" + L("按序列 / 窗口的 LAI", "LAI per sequence / window") + "</summary>" +
                         table(lh, lr, limit=200) + "</details>")

    # Telomeres
    ts = kv_file(AQ("telomere", "telomere_summary.tsv"))
    th, tr = read_tsv(AQ("telomere", "telomere_per_sequence.tsv"))
    if ts or th:
        parts.append("<h3>" + L("端粒（tidk）", "Telomeres (tidk)") + "</h3>" +
                     table([L("项目", "Item"), L("值", "Value")], ts) + table(th, tr, limit=60))

    # FCS-GX
    fh_, fr_ = read_tsv(AQ("fcs_gx", "fcs_gx_report.txt"))
    if fh_:
        parts.append("<h3>" + L("污染筛查（FCS-GX）", "Contamination screen (FCS-GX)") + "</h3><p>" +
                     L("只报告，流程不会自动删除任何序列。", "Reported only; the workflow does not remove any sequence.") +
                     "</p>" + (table(fh_, fr_) if fr_ else "<p>" + L("未发现污染。", "No contamination found.") + "</p>"))

    S.append(("assembly", L("00_Data 基因组组装质量", "00_Data Genome assembly quality"), "".join(parts)))

    # 00_Data: RNA-seq raw data
    fastp = []
    for f_ in sorted(glob.glob(O("00_Data/qc/fastp/*.json"))):
        try:
            import json
            d = json.load(open(f_))
            bf, af = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
            fastp.append([os.path.basename(f_)[:-5], bf["total_reads"], af["total_reads"],
                          f"{100 * af['q30_rate']:.1f}", f"{100 * af['gc_content']:.1f}"])
        except Exception:
            continue
    nano = []
    for f_ in sorted(glob.glob(O("00_Data/qc/nanoplot/*/NanoStats.txt"))):
        kv = {}
        for line in open(f_):
            if ":" in line:
                k, v = line.split(":", 1)
                kv[k.strip()] = v.strip().split()[0] if v.strip() else ""
        nano.append([os.path.basename(os.path.dirname(f_)), kv.get("Number of reads", ""),
                     kv.get("Mean read length", ""), kv.get("Read length N50", ""), kv.get("Mean read quality", "")])
    S.append(("rnaqc", L("00_Data RNA-seq 原始数据", "00_Data RNA-seq raw data"),
              "<h3>" + L("短读（fastp）", "Short reads (fastp)") + "</h3>" +
              table(["sample", "reads_before", "reads_after", "Q30_%", "GC_%"], fastp) +
              "<h3>" + L("长读（NanoPlot）", "Long reads (NanoPlot)") + "</h3>" +
              table(["sample", "reads", "mean_length", "N50", "mean_quality"], nano)))
    # 01_Repeat
    R = lambda n: O("01_Repeat/stats", n + ".tsv")
    rparts = ["<p>" + L("EDTA 构建物种特异 TE 库并做全基因组 TE 注释；RepeatMasker 用该库对基因组软屏蔽"
                        "（同时屏蔽简单重复和低复杂度序列）。重叠的注释只计一次。",
                        "EDTA builds a species-specific TE library and annotates TEs genome-wide; RepeatMasker "
                        "soft-masks the genome with this library (also masking simple repeats and low-complexity "
                        "sequence). Overlapping annotations are counted once.") + "</p>"]
    h, r = read_tsv(R("overview"))
    if h:
        rparts.append("<h3>" + L("总览", "Overview") + "</h3>" + table(h, r))
    h, r = read_tsv(R("te_classes"))
    if h and r:
        tops = [x for x in r if x[1] == "(all)"]
        fig = bar_figure([x[0] for x in tops], [float(x[4]) if x[4] != "-" else 0 for x in tops], "% of genome",
                         figdir, "te_classes", (
                             L("各类转座子占基因组的比例。", "Share of the genome per TE class."),
                             L("EDTA 注释中各 TE 大类（各超家族合计）覆盖的碱基占基因组长度的百分比；重叠部分只计一次。",
                               "Bases covered by each TE class of the EDTA annotation (all superfamilies together) as "
                               "% of the genome length; overlaps are counted once.")),
                         values_text="{:.2f}%")
        rparts.append("<h3>" + L("TE 类别组成", "TE class composition") + "</h3>" + fig + table(h, r))
    h, r = read_tsv(R("te_methods"))
    if h and r:
        rparts.append("<h3>" + L("结构鉴定（完整元件）与同源注释", "Structural identification (intact elements) and "
                                                         "homology annotation") + "</h3><p class='muted'>" +
                      L("structural = 由结构特征鉴定的完整元件；homology = 根据 TE 库同源比对注释的片段。",
                        "structural = intact elements identified by structural features; homology = fragments "
                        "annotated by similarity to the TE library.") + "</p>" + table(h, r))
    h, r = read_tsv(R("library"))
    if h and r:
        n_lib = sum(int(x[1]) for x in r)
        rparts.append("<h3>" + L("TE 库", "TE library") + "</h3><p class='muted'>" +
                      L(f"共 {n_lib} 条一致序列。", f"{n_lib} consensus sequences in total.") + "</p>" + table(h, r))
    h, r = read_tsv(R("landscape"))
    if h and r and plt is not None:
        classes = sorted({x[1] for x in r})
        divs = sorted({int(x[0]) for x in r})
        val = {(int(x[0]), x[1]): float(x[3]) if x[3] != "-" else 0 for x in r}
        f, ax = plt.subplots(figsize=(8, 3.5))
        bottom = [0.0] * len(divs)
        palette = ["#2f6fb0", "#e0a43a", "#5a9e6f", "#c8553d", "#8c6bb1", "#7fb0e6", "#999999", "#d4a5a5"]
        for i, cl in enumerate(classes):
            ys = [val.get((d, cl), 0) for d in divs]
            ax.bar(divs, ys, bottom=bottom, width=1.0, color=palette[i % len(palette)], label=cl)
            bottom = [b + y for b, y in zip(bottom, ys)]
        ax.set_xlabel("Kimura substitution level (%)")
        ax.set_ylabel("% of genome")
        legend_below(ax, ncol=4)                       # outside the axes: the bars fill the plot area
        rparts.append("<h3>" + L("重复序列景观（Kimura 距离）", "Repeat landscape (Kimura distance)") +
                      "</h3><p class='muted'>" + L("距离越小，插入越年轻。", "The smaller the distance, the younger "
                                                                         "the insertion.") + "</p>" +
                      fig_to_img(f, figdir, "repeat_landscape", (
                          L("重复序列景观。", "Repeat landscape."),
                          L("横轴为拷贝与其一致序列之间的 Kimura 距离（经 CpG 校正，%），纵轴为该距离区间内各 TE 类别占基因组的百分比。"
                            "距离越小代表插入越年轻；0 附近的峰常来自近期活跃或刚扩增的家族。",
                            "x axis: Kimura distance between copies and their consensus (CpG-adjusted, %); y axis: "
                            "share of the genome per TE class in each distance bin. Smaller distances mean younger "
                            "insertions; a peak near 0 usually comes from recently active or amplified families."))))
    age_f = R("ltr_age")
    if os.path.exists(age_f):
        lines = [l.rstrip("\n").split("\t") for l in open(age_f) if l.strip()]
        note = [l[0].lstrip("# ") for l in lines if l[0].startswith("#")]
        body = [l for l in lines[1:] if not l[0].startswith("#")][1:] if len(lines) > 2 else []
        title = "<h3>" + L("完整 LTR 插入时间", "Insertion time of intact LTRs") + "</h3>"
        if body:
            by_sf = {}
            for x in body:
                by_sf.setdefault(x[1], []).append(float(x[3]))
            summ = [[sf, len(v), f"{sorted(v)[len(v) // 2]:.2f}", f"{min(v):.2f}", f"{max(v):.2f}"]
                    for sf, v in sorted(by_sf.items(), key=lambda t: -len(t[1]))]
            fig = ""
            n_ltr = sum(len(v) for v in by_sf.values())
            if plt is not None and n_ltr < 20:
                fig = "<p class='muted'>" + L(f"只有 {n_ltr} 个完整 LTR，数量太少，不画分布图，见下表。",
                                              f"Only {n_ltr} intact LTRs, too few for a distribution plot; see the "
                                              "table.") + "</p>"
            elif plt is not None:
                f, ax = plt.subplots(figsize=(7, 3.2))
                for i, (sf, v) in enumerate(sorted(by_sf.items(), key=lambda t: -len(t[1]))[:4]):
                    ax.hist(v, bins=40, alpha=0.6, label=sf, color=["#2f6fb0", "#e0a43a", "#5a9e6f", "#c8553d"][i])
                ax.set_xlabel("insertion time (Mya)")
                ax.set_ylabel("intact LTRs")
                ax.legend(frameon=False)
                fig = fig_to_img(f, figdir, "ltr_insertion_time", (
                    L("完整 LTR 反转录转座子的插入时间。", "Insertion time of intact LTR retrotransposons."),
                    L("按超家族（最多 4 个）分别统计；插入时间 T = K / (2u)，K 为两端 LTR 之间的分化度，"
                      "u 为替换速率（EDTA -u，见配置）。",
                      "Per superfamily (at most 4); insertion time T = K / (2u), K = divergence between the two LTRs, "
                      "u = substitution rate (EDTA -u, see the configuration).")))
            rparts.append(title + "<p class='muted'>" + html.escape("; ".join(note)) + "</p>" + fig +
                          table([L("超家族", "Superfamily"), L("完整 LTR 数", "Intact LTRs"),
                                 L("中位数 (Mya)", "Median (Mya)"), L("最小", "Min"), L("最大", "Max")], summ))
        else:
            rparts.append(title + "<p class='muted'>" + L("未找到完整 LTR 列表或其中没有插入时间。",
                                                          "No intact LTR list found, or it has no insertion times.") +
                          "</p>")
    sumf = glob.glob(O("01_Repeat/EDTA/*.TEanno.sum"))
    if sumf:
        rparts.append("<details><summary>" + L("EDTA 注释汇总原文", "EDTA annotation summary (original text)") +
                      "</summary><pre>" + html.escape(open(sumf[0]).read()[:20000]) + "</pre></details>")
    S.append(("repeat", L("01_Repeat 重复序列", "01_Repeat Repeats"), "".join(rparts)))
    # 02_RNA: alignment
    ah, ar = read_tsv(O("02_RNA/qc/alignment_summary.tsv"))
    mh, mr = read_tsv(O("02_RNA/qc/merged_bam_summary.tsv"))
    fig = ""
    if plt is not None and ah and ar:
        i_m, i_s = ah.index("mapping_%"), ah.index("spliced_%")
        labels = [f"{x[0]} ({x[2]})" for x in ar]
        f, ax = plt.subplots(figsize=(7.5, 0.45 * len(ar) + 1.2))
        y = list(range(len(ar)))
        ax.barh([v - 0.2 for v in y], [float(x[i_m]) for x in ar], height=0.4, color="#2f6fb0", label="mapped %")
        ax.barh([v + 0.2 for v in y], [float(x[i_s]) for x in ar], height=0.4, color="#e0a43a", label="spliced % of mapped")
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_xlim(0, 100)
        ax.set_xlabel("%")
        legend_below(ax, ncol=2)
        fig = fig_to_img(f, figdir, "rna_alignment", (
            L("RNA-seq 比对率。", "RNA-seq alignment rates."),
            L("每行一个样本（括号内为平台）。蓝色 = 比对上的 reads 占总 reads 的百分比；"
              "橙色 = 比对上的 reads 中含剪接（CIGAR 中有 N）的比例，长读通常远高于短读。",
              "One row per sample (platform in brackets). Blue = aligned reads as % of all reads; orange = share of "
              "aligned reads that are spliced (N in the CIGAR), usually much higher for long reads than short reads.")))
    S.append(("rnaalign", L("02_RNA RNA-seq 比对", "02_RNA RNA-seq alignment"),
              "<p>" + L("短读用 HISAT2 比对，长读用 minimap2 比对（预设按平台选择）。比对率 = 比对上的 reads / 总 reads"
                        "（samtools stats，不含 secondary/supplementary）；unique/multi 来自 HISAT2 汇总；"
                        "spliced % = 含内含子（CIGAR 中有 N）的主比对占比对上 reads 的比例；error rate 为错配碱基 / 比对碱基。",
                        "Short reads are aligned with HISAT2, long reads with minimap2 (preset by platform). Alignment "
                        "rate = aligned reads / all reads (samtools stats, without secondary/supplementary); "
                        "unique/multi from the HISAT2 summary; spliced % = primary alignments with an intron (N in "
                        "the CIGAR) as share of aligned reads; error rate = mismatched bases / aligned bases.") + "</p>"
              "<h3>" + L("各样本", "Samples") + "</h3>" + fig + table(ah, ar) +
              "<h3>" + L("合并 BAM 及其下游用途", "Merged BAMs and their use downstream") + "</h3>" + table(mh, mr) +
              long_read_weights_html(O("02_RNA/qc/long_read_weights.tsv"))))
    # 02-04: annotation sources (same layout for EviAnn, BRAKER4, GeMoMa)
    prot = EVIANN_PROTEINS.get(a.eviann_proteins, (a.eviann_proteins, a.eviann_proteins))
    prot = prot[0] if LANG == "zh" else prot[1]
    S.append(("eviann", L("02_RNA EviAnn 注释", "02_RNA EviAnn annotation"), annotation_section(
        O, "eviann", O("02_RNA/EviAnn/stats/eviann_stats.tsv"), O("02_RNA/EviAnn/stats/eviann_distributions.tsv"),
        O("02_RNA/EviAnn/qc"), figdir,
        "<p>" + L(f"EviAnn 以 RNA 比对为证据进行注释，每个样本单独一行输入（<b>{html.escape(a.eviann_tag)}</b>），"
                  f"EviAnn 对每个样本分别组装再合并；确定 CDS 所用的蛋白集：<b>{html.escape(prot)}</b>。"
                  "结果按是否含 CDS 拆分为编码基因（进入 EVM）和 lncRNA（单独输出；EviAnn 只保留至少在 2 个样本中检测到的 lncRNA）。",
                  f"EviAnn annotates from the RNA alignments, one input line per sample (<b>{html.escape(a.eviann_tag)}"
                  "</b>), assembling each sample separately and then merging; protein set used to call CDSs: "
                  f"<b>{html.escape(prot)}</b>. The result is split by CDS into coding genes (to EVM) and lncRNAs "
                  "(written separately; EviAnn keeps a lncRNA only when it is found in at least 2 samples).") + "</p>")))

    sp_h, sp_r = read_tsv(O("03_Ab_initio/stats/braker4_support.tsv"))
    sp_note = ""
    sp_file = O("03_Ab_initio/stats/braker4_support.tsv")
    if os.path.exists(sp_file):
        sp_note = "".join(f"<p class='muted'>{html.escape(l[1:].strip())}</p>" for l in open(sp_file) if l.startswith("#"))
    sp_r = [r for r in sp_r if not r[0].startswith("#")]
    mode = html.escape(a.braker4_mode)
    S.append(("braker4", L("03_Ab_initio BRAKER4 注释", "03_Ab_initio BRAKER4 annotation"), annotation_section(
        O, "braker4", O("03_Ab_initio/stats/braker4_stats.tsv"), O("03_Ab_initio/stats/braker4_distributions.tsv"),
        O("03_Ab_initio/qc"), figdir,
        "<div class='note'>" +
        L(f"BRAKER4 模式：<b>{mode}</b>（按数据自动选择：短读 + PacBio = dual，有短读无 PacBio = ETP，仅 PacBio = IsoSeq，"
          "仅 ONT 或无 RNA = EP；ONT 长读不进入 BRAKER4）。输入为软屏蔽基因组、OrthoDB 与同源物种蛋白。",
          f"BRAKER4 mode: <b>{mode}</b> (chosen from the data: short + PacBio = dual, short without PacBio = ETP, "
          "PacBio only = IsoSeq, ONT only or no RNA = EP; ONT long reads are not passed to BRAKER4). Inputs: the "
          "soft-masked genome, OrthoDB and the homologous species' proteins.") + "</div>",
        extra="<h3>" + L("证据支持（gene_support.tsv）", "Evidence support (gene_support.tsv)") + "</h3>" + sp_note +
              "<p class='muted'>" + L("BRAKER4 未说明该文件的列含义，此处逐列统计非零/为真的模型比例。",
                                      "BRAKER4 does not document the columns of this file; for each column the share "
                                      "of models with a non-zero / true value is given.") + "</p>" +
              table(sp_h, sp_r))))

    rf_h, rf_r = read_tsv(O("04_Homology/stats/gemoma_references.tsv"))
    at_h, at_r = read_tsv(O("04_Homology/stats/gemoma_attributes.tsv"))
    mp_h, mp_r = read_tsv(O("04_Homology/stats/miniprot_species.tsv"))
    sp_qc = []
    for h_ in a.homologs:
        for k, v in protein_qc_rows(O("04_Homology/qc", f"GeMoMa_{h_}")):
            sp_qc.append([h_, k, v])
    S.append(("gemoma", L("04_Homology 同源注释（GeMoMa / miniprot）", "04_Homology Homology annotation (GeMoMa / "
                                                                    "miniprot)"), annotation_section(
        O, "gemoma", O("04_Homology/stats/gemoma_stats.tsv"), O("04_Homology/stats/gemoma_distributions.tsv"),
        O("04_Homology/qc/GeMoMa_merged"), figdir,
        "<p>" + L("每个同源物种单独运行一次 GeMoMa：把该物种的基因模型投射到目标基因组，并用 RNA-seq 的剪接位点和覆盖度"
                  "（ERE，只提取一次）校正。各物种的预测再用 GeMoMa GAF + AnnotationFinalizer 合并为 GeMoMa_merged，"
                  "作为 EVM 的输入。下表中 GeMoMa_&lt;物种&gt; 为各物种单独的结果；蛋白集质量一栏为合并结果。",
                  "GeMoMa runs once per homologous species: the species' gene models are projected onto the target "
                  "genome and corrected with RNA-seq splice sites and coverage (ERE, extracted once). The predictions "
                  "of all species are combined with GeMoMa GAF + AnnotationFinalizer into GeMoMa_merged, the EVM "
                  "input. GeMoMa_&lt;species&gt; in the table is the result of one species; the protein set quality "
                  "refers to the merged result.") + "</p>",
        extra="<h3>" + L("各物种单独结果的蛋白集质量", "Protein set quality per species") + "</h3>" +
              (table([L("物种", "Species"), L("评估", "Assessment"), L("结果", "Result")], sp_qc) if sp_qc else
               "<p class='muted'>" + L("未做按物种的评估（qc.per_homolog: false）。",
                                       "No per-species assessment (qc.per_homolog: false).") + "</p>") +
              "<h3>" + L("各物种对合并结果的贡献", "Contribution of each species to the merged result") +
              "</h3><p class='muted'>" +
              L("merged_representative = 合并结果中代表预测来自该物种的转录本；"
                "merged_supported = 该物种也预测出相同 CDS 的转录本（含 GAF alternative）。",
                "merged_representative = transcripts of the merged result whose representative prediction comes from "
                "the species; merged_supported = transcripts for which the species also predicted the same CDS "
                "(including GAF alternatives).") + "</p>" + table(rf_h, rf_r) +
              "<h3>" + L("合并结果的证据属性", "Evidence attributes of the merged result") + "</h3><p class='muted'>" +
              L("evidence = 预测出该转录本的物种数；tie = 有 RNA 支持的内含子比例；"
                "tpc = 被 RNA 覆盖的转录本比例（后两项只在使用 RNA 证据时有意义）。",
                "evidence = number of species predicting the transcript; tie = share of introns supported by RNA; "
                "tpc = share of the transcript covered by RNA (the last two only mean something with RNA evidence).") +
              "</p>" + table(at_h, at_r) +
              "<h3>" + L("miniprot 蛋白比对（按物种）", "miniprot protein alignments (per species)") +
              "</h3><p class='muted'>" +
              L("每个物种的蛋白单独比对到基因组；合并后作为 EVM 的 PROTEIN 证据，也用于修正（R 规则）和权重调优。",
                "The proteins of each species are aligned to the genome separately; combined, they are EVM's PROTEIN "
                "evidence and are also used by the refinement (R rules) and the weight tuning.") + "</p>" +
              table(mp_h, mp_r))))

    # QC comparison
    h, r = read_tsv(O("05_EVM/qc/summary.tsv"))
    ref_html = ""
    ref_f = O("05_EVM/qc/reference_comparison.tsv")
    if os.path.exists(ref_f):
        lines = [l.rstrip("\n").split("\t") for l in open(ref_f) if l.strip()]
        note = [l[0].lstrip("# ") for l in lines if l[0].startswith("#")]
        body = [l for l in lines if not l[0].startswith("#")]
        if body:
            rh, rr = body[0], body[1:]
            fig = ""
            if plt is not None and rr:
                f, ax = plt.subplots(figsize=(8, 0.5 * len(rr) + 1.2))
                y = list(range(len(rr)))
                for k, (col, colr) in enumerate((("gene_recall_%", "#2f6fb0"), ("intron_recall_%", "#5a9e6f"),
                                                  ("tx_recall_%", "#e0a43a"))):
                    vals = [float(x[rh.index(col)]) if x[rh.index(col)] != "-" else 0 for x in rr]
                    ax.barh([v + (k - 1) * 0.27 for v in y], vals, height=0.27, color=colr, label=col.replace("_%", ""))
                ax.set_yticks(y)
                ax.set_yticklabels([x[0] for x in rr])
                ax.invert_yaxis()
                ax.set_xlim(0, 100)
                ax.set_xlabel("% of reference")
                legend_below(ax, ncol=3)
                fig = fig_to_img(f, figdir, "reference_comparison", (
                    L("与参考注释的一致性（召回）。", "Agreement with the reference annotation (recall)."),
                    L("gene_recall = 参考基因中与该注释集有同链 CDS 重叠的比例；intron_recall = 参考的 CDS 内含子"
                      "被完全相同地预测出的比例；tx_recall = 参考转录本的完整 CDS 结构被完全相同地预测出的比例。精确度见下表。"
                      "参考注释本身也是预测结果时，反映的是一致程度而非准确率。",
                      "gene_recall = reference genes with a CDS overlap on the same strand; intron_recall = reference "
                      "CDS introns predicted exactly; tx_recall = reference transcripts whose complete CDS structure is "
                      "predicted exactly. Precision is in the table. When the reference is itself a prediction, these "
                      "values measure agreement, not accuracy.")))
            ref_html = ("<h3>" + L("与参考注释的一致性（CDS 水平）", "Agreement with the reference annotation (CDS "
                                                                    "level)") + "</h3>" + ref_compare_help() +
                        "<p class='muted'>" + html.escape("; ".join(note)) + "</p>" + fig + table(rh, rr))
    ch, cr = read_tsv(O("05_EVM/stats/all_sets_stats.tsv"))
    S.append(("qc", L("各注释来源对比（BUSCO / OMArk / 基因结构）", "Comparison of annotation sources (BUSCO / OMArk / "
                                                               "gene structure)"),
              "<p>" + L("所有来源使用同一组用户指定的 BUSCO lineage，BUSCO / OMArk 都用每个基因最长的蛋白。各行：",
                        "All sources use the same user-specified BUSCO lineages; BUSCO / OMArk use the longest protein "
                        "of every gene. Rows:") + "</p><ul>" +
              L("<li><b>EviAnn / BRAKER4 / GeMoMa_*</b>：三路预测（RNA 注释、从头预测、同源投射），即 EVM 的输入；"
                "GeMoMa_&lt;物种&gt; 为单个同源物种，GeMoMa_merged 为合并结果</li>"
                "<li><b>EVM</b>：EVM 整合的共识模型（没有 UTR，每个基因一个转录本）</li>"
                "<li><b>refine</b>：与 EviAnn 比较修正之后（引入 UTR 和异构体，补入、替换基因），过滤之前，含全部异构体</li>"
                "<li><b>final</b>：TE / 长度过滤和重新命名之后的最终结果；主结果每个基因一个转录本。"
                "与 refine 在 BUSCO / OMArk 上的差别只来自过滤删掉的基因</li>",
                "<li><b>EviAnn / BRAKER4 / GeMoMa_*</b>: the three prediction routes (RNA-based, ab initio, homology "
                "projection), the inputs of EVM; GeMoMa_&lt;species&gt; = one homologous species, GeMoMa_merged = "
                "combined</li>"
                "<li><b>EVM</b>: EVM consensus models (no UTRs, one transcript per gene)</li>"
                "<li><b>refine</b>: after the refinement against EviAnn (UTRs and isoforms added, genes added or "
                "replaced), before filtering, with all isoforms</li>"
                "<li><b>final</b>: the final result after the TE / length filter and renaming; the main result has one "
                "transcript per gene. Its BUSCO / OMArk differ from refine only by the genes removed by the filter</li>") +
              "</ul>" +
              busco_figure(h, r, figdir, "busco_comparison", (
                  L("各注释来源的蛋白集 BUSCO。", "Protein-set BUSCO of the annotation sources."),
                  L("每行一个注释来源 × 一个 BUSCO 数据集，输入为各来源每个基因的最长蛋白；"
                    "按 S / D / F / M 堆叠，右侧为完整度 C。各来源使用同一组数据集，可以直接比较。",
                    "One row per annotation source x BUSCO dataset, using the longest protein of every gene; stacked "
                    "S / D / F / M with completeness C on the right. All sources use the same datasets and can be "
                    "compared directly."))) + omark_figure(h, r, figdir) + table(h, r) +
              "<h3>" + L("基因结构对比", "Gene structure comparison") + "</h3><p class='muted'>" +
              L("EviAnn 一列只含编码基因；final 为过滤和重新命名后的最终基因集。",
                "The EviAnn column contains coding genes only; final is the final gene set after filtering and "
                "renaming.") + "</p>" + table(ch, cr) + ref_html))
    # 05_EVM: EVM integration
    impl = kv_file(O("05_EVM/logs/evm_impl_check.txt"))
    impl_txt = ""
    if impl:
        bits = [impl[0][0]] + impl[0][1].split("\t")
        impl_txt = ("<div class='note'>" + L("EVM 实现：", "EVM implementation: ") +
                    f"<b>{html.escape(bits[0])}</b>" + L(f"（{html.escape(bits[-1])}）", f" ({html.escape(bits[-1])})") +
                    f"<br><span class='muted'>{html.escape(bits[1] if len(bits) > 2 else '')}</span></div>")
    w_rows = [l.rstrip("\n").split("\t") for l in open(O("05_EVM/inputs/weights.txt"))
              if l.strip()] if os.path.exists(O("05_EVM/inputs/weights.txt")) else []
    S.append(("evm", L("05_EVM EVM 整合", "05_EVM EVM integration"), annotation_section(
        O, "evm", O("05_EVM/stats/evm_stats.tsv"), O("05_EVM/stats/evm_distributions.tsv"), O("05_EVM/qc/evm"),
        figdir, impl_txt + tune_html(O("05_EVM/tune/tuning_results.tsv"), O("05_EVM/tune/best_weights.txt"), figdir)
        + "<h3>" + L("输入权重", "Input weights") + "</h3>" +
        table([L("类别", "Class"), L("来源", "Source"), L("权重", "Weight")], w_rows),
        extra=source_support_html(O("05_EVM/stats/evm_source_support.tsv"),
                                  O("05_EVM/stats/evm_source_combinations.tsv"), figdir, "evm"))))

    # refine
    dec_f = O("05_EVM/refine/refine_decisions.tsv")
    notes, h, r = [], None, []
    if os.path.exists(dec_f):
        lines = [l.rstrip("\n").split("\t") for l in open(dec_f) if l.strip()]
        notes = [l[0].lstrip("# ") for l in lines if l[0].startswith("#")]
        body = [l for l in lines if not l[0].startswith("#")]
        if body:
            h, r = body[0], body[1:]
    rparts = [L("<p>EVM 与 EviAnn 编码基因按 CDS 在同一条链上的重叠划分为位点，按 CDS 内含子链比较：</p><ul>"
                "<li>R1 一致：保留 EVM 的 CDS，引入 EviAnn 的 UTR 和其它完整异构体</li>"
                "<li>R2 EVM 缺失：EviAnn 基因有完整 ORF、RNA 支持且 CDS 与 TE 重叠不超过阈值时补入</li>"
                "<li>R3 结构冲突：EviAnn 内含子链被 RNA 支持而 EVM 不被支持时替换；两者都被支持时把 EviAnn 模型加为异构体</li>"
                "<li>R4 拆分 / 融合：有 RNA 支持时采用 EviAnn 的结构，有 reads 连接两个基因时保留融合；证据不足时保留 EVM 并标记 needs_review</li>"
                "<li>RNA 支持：有长读时为长读内含子链与模型完全一致；只有短读时，1 个内含子看剪接 reads（B：绝对数与相对竞争内含子的比例），"
                "≥ 2 个内含子还要求 StringTie 组装出相同的内含子链（A 且 B），单外显子模型看单外显子 StringTie 转录本（C）；"
                "决策表 support_basis 列记录依据</li>"
                "<li>R5 EVM 独有：检查蛋白、RNA、BRAKER4 证据，没有证据的标记 refine_support=none</li>"
                "<li>R6 lncRNA：单独输出，不参与合并</li></ul>",
                "<p>EVM and EviAnn coding genes are grouped into loci by CDS overlap on the same strand and compared "
                "by their CDS intron chains:</p><ul>"
                "<li>R1 identical: EVM's CDS is kept, EviAnn's UTRs and other complete isoforms are added</li>"
                "<li>R2 missing in EVM: an EviAnn gene is added when it has a complete ORF, RNA support and a CDS "
                "overlap with TEs below the threshold</li>"
                "<li>R3 structural conflict: replaced when EviAnn's intron chain is supported by RNA and EVM's is not; "
                "when both are supported the EviAnn model is added as an isoform</li>"
                "<li>R4 split / fusion: EviAnn's structure is taken with RNA support; a fusion is kept when reads "
                "connect the two genes; with insufficient evidence EVM is kept and marked needs_review</li>"
                "<li>RNA support: with long reads, a long-read intron chain identical to the model; with short reads "
                "only, junction reads for 1 intron (B: absolute count and ratio to competing introns), plus an "
                "identical StringTie intron chain for ≥ 2 introns (A and B), and a single-exon StringTie transcript for "
                "single-exon models (C); the support_basis column of the decision table records the basis</li>"
                "<li>R5 EVM only: protein, RNA and BRAKER4 evidence is checked; genes without evidence are marked "
                "refine_support=none</li>"
                "<li>R6 lncRNA: written separately, not merged</li></ul>") + refine_rules_html(a.config)]
    rparts += [f"<p class='muted'>{html.escape(n_)}</p>" for n_ in notes]
    if h:
        c = Counter((x[h.index("rule")], x[h.index("action")]) for x in r)
        keys = sorted(c, key=lambda k: (k[0], -c[k]))
        cg = Counter()                                   # merge_utr+N_isoforms -> one bar
        for (ru, ac), n_ in c.items():
            cg[(ru, re.sub(r"\+\d+_isoforms$", " (+isoforms)", ac))] += n_
        gkeys = sorted(cg, key=lambda k: (k[0], -cg[k]))
        rule_col = {"R1": "#2f6fb0", "R2": "#5a9e6f", "R3": "#e0a43a", "R5": "#8c6bb1"}
        cols = [rule_col.get(k[0][:2], "#c8553d") for k in gkeys]
        rparts.append(bar_figure(
            [f"{k[0]} {k[1]}" for k in gkeys], [cg[k] for k in gkeys], "loci", figdir, "refine_rules",
            (L("修正规则的处理结果。", "Outcome of the refinement rules."),
             L("每个位点（EVM 与 EviAnn 编码基因按同链 CDS 重叠分组）按适用的规则和处理方式计数；"
               "merge_utr (+isoforms) 合并了引入不同数目异构体的情况，明细见下表。颜色区分规则。",
               "Loci (EVM and EviAnn coding genes grouped by CDS overlap on the same strand) counted by the rule "
               "applied and the action taken; merge_utr (+isoforms) combines the cases with different numbers of "
               "added isoforms (details in the table). Colours distinguish the rules.")),
            colors=cols, legend=[(k, v) for k, v in (("R1", "#2f6fb0"), ("R2", "#5a9e6f"), ("R3", "#e0a43a"),
                                                       ("R4", "#c8553d"), ("R5", "#8c6bb1"))
                                 if any(x[0].startswith(k) for x in gkeys)], values_text="{:.0f}"))
        rparts.append(table([L("规则", "Rule"), L("处理", "Action"), L("位点数", "Loci")],
                            [[k[0], k[1], c[k]] for k in keys]))
        review = [x for x in r if "needs_review" in x[h.index("note")]]
        if review:
            rparts.append("<h3>" + L(f"需要人工检查的位点（{len(review)} 个）", f"Loci to check manually ({len(review)})") +
                          "</h3>" + table(h, review, limit=50))
        rparts.append("<details><summary>" + L("全部决策（前 200 行）", "All decisions (first 200 rows)") +
                      "</summary>" + table(h, r, limit=200) + "</details>")
    S.append(("refine", L("05_EVM 修正（EVM vs EviAnn）", "05_EVM Refinement (EVM vs EviAnn)"), "".join(rparts)))

    # filtering
    h, r = read_tsv(O("05_EVM/final/removed_genes.tsv"))
    fparts = []
    if h:
        c = Counter(x[1].split("_")[0] for x in r)
        fparts.append(table([L("原因", "Reason"), L("基因数", "Genes")], sorted(c.items())))
    dec = O("05_EVM/final/filter_decisions.tsv")
    if os.path.exists(dec):
        rows_ = [l.rstrip("\n").split("\t") for l in open(dec) if l.strip()]
        notes = [x[0].lstrip("# ") for x in rows_ if x[0].startswith("#")]
        body = [x for x in rows_ if not x[0].startswith("#")]
        if body:
            hd, body = body[0], body[1:]
            ix = {k: i for i, k in enumerate(hd)}
            thr = None
            for n_ in notes:
                m = re.search(r"max_te_overlap=([\d.]+)", n_)
                if m:
                    thr = float(m.group(1))
            high_te = [x for x in body if thr is not None and float(x[ix["te_overlap_frac"]]) > thr]
            rescued = [x for x in high_te if x[ix["rna_supported"]] == "yes"]
            by = Counter()
            for x in rescued:
                for k in ("rna_stringtie", "rna_eviann", "rna_braker4"):
                    if x[ix[k]] == "yes":
                        by[k] += 1
            summ = [[L("参与过滤的基因", "Genes checked by the filter"), len(body)],
                    [L(f"TE 重叠 > {thr}", f"TE overlap > {thr}"), len(high_te)],
                    [L("其中因有 RNA 证据而保留", "of these kept because of RNA evidence"), len(rescued)],
                    [L("　StringTie 支持", "  StringTie support"), by["rna_stringtie"]],
                    [L("　EviAnn 支持", "  EviAnn support"), by["rna_eviann"]],
                    [L("　BRAKER4 gene_support 支持", "  BRAKER4 gene_support support"), by["rna_braker4"]],
                    [L("有 TEsorter 命中（仅记录）", "TEsorter hits (recorded only)"),
                     sum(x[ix["tesorter"]] != "-" for x in body)],
                    [L("CDS 短于下限但有证据而保留（short_cds=yes）", "CDS below the minimum but kept for evidence "
                                                                  "(short_cds=yes)"),
                     sum(x[ix["reason"]].startswith("short_CDS") for x in body)]]
            fparts.append("<h3>" + L("TE 过滤明细", "TE filter details") + "</h3>" +
                          "".join(f"<p class='muted'>{html.escape(n_)}</p>" for n_ in notes) +
                          table([L("项目", "Item"), L("数目", "Count")], summ))
    S.append(("filter", L("05_EVM 过滤", "05_EVM Filtering"),
              "<p>" + L("基因 CDS 与 TE 重叠超过阈值且没有 RNA 证据时才剔除；CDS 短于 refine.min_cds_aa 的模型单独剔除，"
                        "但有 RNA 证据（无 RNA-seq 时为蛋白证据）的保留并标记 short_cds=yes（refine.keep_short_supported）；"
                        "TEsorter 只记录。",
                        "A gene is removed only when its CDS overlap with TEs exceeds the threshold and it has no RNA "
                        "evidence; models with a CDS below refine.min_cds_aa are removed separately, except those with RNA "
                        "evidence (without RNA-seq: protein evidence), which are kept and marked short_cds=yes "
                        "(refine.keep_short_supported); TEsorter hits are only recorded.") + "</p>" + "".join(fparts)))

    # final gene set
    ch, cr = read_tsv(O("05_EVM/stats/all_sets_stats.tsv"))
    prog = []
    if ch:
        idx = {m[0]: m for m in cr}
        for col in ("EVM", "refine", "final"):
            if col in ch:
                j = ch.index(col)
                prog.append([col, idx.get("genes", ["", ""] * 9)[j], idx.get("transcripts", ["", ""] * 9)[j]])
    files = []
    for f_ in sorted(glob.glob(O("05_EVM/final", a.prefix + ".*"))):
        size = os.path.getsize(f_)
        if f_.endswith((".fa", ".faa", ".fasta")):
            n_rec = sum(1 for l in open(f_) if l.startswith(">"))
            kind = L("序列", "sequences")
        else:
            n_rec = sum(1 for l in open(f_) if "\tgene\t" in l)
            kind = L("基因", "genes")
        files.append([os.path.relpath(f_, a.outdir), f"{size / 1e6:.2f} MB", f"{n_rec} {kind}"])
    pfx = html.escape(a.prefix)
    S.append(("finalset", L("05_EVM 最终注释结果", "05_EVM Final annotation"), annotation_section(
        O, "final", O("05_EVM/stats/final_stats.tsv"), O("05_EVM/stats/final_distributions.tsv"),
        O("05_EVM/qc/final"), figdir,
        "<p>" + L("最终基因集：EVM 整合 → 与 EviAnn 比较修正 → TE / 长度过滤 → 重新命名并还原原始序列名。"
                  "lncRNA 来自 EviAnn，单独输出。",
                  "Final gene set: EVM integration → refinement against EviAnn → TE / length filter → renaming and "
                  "restoring the original sequence names. lncRNAs come from EviAnn and are written separately.") +
        "</p><div class='note'>" +
        L(f"主结果（<code>{pfx}.gff3</code> 及 cds/pep/cdna）每个基因只保留一个转录本："
          "CDS 最长者，CDS 等长时取转录本最长者；本节统计、BUSCO/OMArk 和功能注释都基于主结果。"
          f"包含全部异构体的结果另存为 <code>{pfx}.all_isoforms.*</code>。",
          f"The main result (<code>{pfx}.gff3</code> and cds/pep/cdna) keeps one transcript per gene: the one with the "
          "longest CDS, with equal CDS the longest transcript; the statistics of this section, BUSCO/OMArk and the "
          f"functional annotation are based on the main result. The result with all isoforms is in "
          f"<code>{pfx}.all_isoforms.*</code>.") + "</div>" +
        (table([L("项目", "Item"), L("数目", "Count")],
               [[k, v] for k, v in kv_file(O("05_EVM/final/isoform_selection.tsv")) if k != "item"] +
               ([["transcripts_with_internal_stop (kept, internal_stop=yes)",
                  str(max(0, len(open(O("05_EVM/final/internal_stops.tsv")).read().splitlines()) - 1))]]
                if os.path.exists(O("05_EVM/final/internal_stops.tsv")) else []))
         if os.path.exists(O("05_EVM/final/isoform_selection.tsv")) else "") +
        "<h3>" + L("各阶段基因数", "Genes per stage") + "</h3>" +
        table([L("阶段", "Stage"), L("基因", "Genes"), L("转录本", "Transcripts")], prog),
        extra=source_support_html(O("05_EVM/stats/final_source_support.tsv"),
                                  O("05_EVM/stats/final_source_combinations.tsv"), figdir, "final") +
              naming_html(O("05_EVM/final/id_map.tsv"), a.prefix) +
              "<h3>" + L("输出文件", "Output files") + "</h3>" +
              table([L("文件", "File"), L("大小", "Size"), L("记录数", "Records")], files))))

    # function
    F = lambda n: O("06_Function/stats", n + ".tsv")
    fparts = ["<p>" + L("对最终基因集主结果的蛋白（每个基因一个转录本）进行 InterProScan 和 eggNOG-mapper 注释，"
                        "统计均以基因为单位。完整的 InterProScan 结果见 06_Function/interproscan/。",
                        "The proteins of the main result of the final gene set (one transcript per gene) are annotated "
                        "with InterProScan and eggNOG-mapper; all statistics count genes. The complete InterProScan "
                        "result is in 06_Function/interproscan/.") + "</p>"]
    oh, orows = read_tsv(F("overview"))
    if oh:
        keep = [x for x in orows if x[0] in ("InterProScan only", "both", "eggNOG-mapper only", "neither (unannotated)")]
        fig = bar_figure([x[0] for x in keep], [float(x[2]) for x in keep], "% genes", figdir, "function_overview", (
            L("功能注释来源。", "Sources of the functional annotation."),
            L("最终基因中被 InterProScan 与 eggNOG-mapper 都注释、只被其一注释、或都没有注释的基因占比。",
              "Share of final genes annotated by both InterProScan and eggNOG-mapper, by only one of them, or by "
              "neither.")),
            values_text="{:.1f}%")
        fparts.append("<h3>" + L("总览", "Overview") + "</h3>" + table(oh, orows) + fig)
    h, r = read_tsv(O("06_Function/merged/function_stats.tsv"))
    if h:
        # Preferred_name / Description are eggNOG name and description fields written to the GFF3
        # (gene_name= / product=), not annotation categories: table only, not in the figure
        fields = {"Preferred_name": L("基因名（gene_name，eggNOG Preferred_name）",
                                      "gene name (gene_name, eggNOG Preferred_name)"),
                  "Description": L("功能描述（product，eggNOG Description）", "description (product, eggNOG Description)")}
        cats = [x for x in r if x[0] not in fields]
        rows_ = [[fields.get(x[0], x[0])] + x[1:] for x in r]
        fparts.append("<h3>" + L("各类注释的基因数", "Genes per annotation category") + "</h3>" +
                      bar_figure([x[0] for x in cats], [float(x[2]) for x in cats], "% genes", figdir, "function_rates", (
                          L("各类功能注释的覆盖率。", "Coverage of the annotation categories."),
                          L("有该类注释的基因占最终基因的百分比。any = 至少有一类；InterPro / Pfam 来自 InterProScan，"
                            "COG、KEGG、EC 来自 eggNOG-mapper，GO 为两者合并（灰色）。基因名和功能描述不是注释类别，只列在下表。",
                            "Genes with the category as % of final genes. any = at least one category; InterPro / Pfam "
                            "from InterProScan, COG, KEGG, EC from eggNOG-mapper, GO from both combined (grey). Gene "
                            "names and descriptions are not categories and are listed in the table only.")),
                          colors=[FUNC_SRC_COL.get(x[0], "#999999") for x in cats],
                          legend=[("InterProScan", "#2f6fb0"), ("eggNOG-mapper", "#e0a43a"), ("both (merged)", "#999999")],
                          values_text="{:.1f}%") +
                      table(h, rows_) +
                      "<p class='muted'>" +
                      L("基因名和功能描述由 eggNOG-mapper 从直系同源基因转移而来，分别写入最终 GFF3 "
                        "（06_Function/merged/&lt;前缀&gt;.func.gff3）的 gene_name= 和 product= 属性。基因名只在同源基因有公认名称时才给出，"
                        "覆盖率低是正常的；功能描述在 eggNOG 找到同源组时一般都有，覆盖率通常与 COG 相同。",
                        "Gene names and descriptions are transferred by eggNOG-mapper from orthologs and written to the "
                        "gene_name= and product= attributes of the final GFF3 (06_Function/merged/&lt;prefix&gt;"
                        ".func.gff3). A gene name is given only when the orthologs have an established name, so low "
                        "coverage is normal; a description is usually present whenever eggNOG finds an orthologous "
                        "group, with coverage similar to COG.") + "</p>")
    h, r = read_tsv(F("interpro_databases"))
    if h:
        fparts.append("<h3>" + L("InterProScan 各成员数据库", "InterProScan member databases") + "</h3>" +
                      bar_figure([x[0] for x in r if float(x[2]) > 0], [float(x[2]) for x in r if float(x[2]) > 0],
                                 "% genes", figdir, "interpro_databases", (
                                     L("InterProScan 各成员数据库的命中率。", "Hit rates of the InterProScan member databases."),
                                     L("被该数据库命中的基因占最终基因的百分比（没有命中的数据库不画）。",
                                       "Genes hit by the database as % of final genes (databases without hits are not "
                                       "shown).")),
                                 values_text="{:.1f}%") +
                      table(h, r))
    for n_, title in (("top_pfam", L("最常见的 Pfam 结构域", "Most frequent Pfam domains")),
                      ("top_interpro", L("最常见的 InterPro 条目", "Most frequent InterPro entries"))):
        h, r = read_tsv(F(n_))
        if h and r:
            fparts.append(f"<h3>{title}</h3>" + table(h, r))
    gh, gr = read_tsv(F("go_namespace"))
    th, tr = read_tsv(F("go_top"))
    go_html = ""
    if gh and gr:
        go_html += table(gh, gr)
    elif th:
        go_html += "<p class='muted'>" + L("未提供 GO 本体文件（function.go_obo），只显示 GO 编号，不分 BP/MF/CC。",
                                           "No GO ontology file (function.go_obo): GO IDs only, without BP/MF/CC.") + "</p>"
    if th and tr:
        go_html += table(th, tr)
    if go_html:
        fparts.append("<h3>GO</h3>" + go_html)
    kh, kr = read_tsv(F("kegg"))
    ph, pr = read_tsv(F("kegg_top_pathways"))
    if kh:
        fparts.append("<h3>KEGG</h3>" + table(kh, kr) +
                      ("<p class='muted'>" + L("最常见的通路（map 编号）：", "Most frequent pathways (map IDs):") + "</p>" +
                       table(ph, pr) if pr else ""))
    ch_, cr_ = read_tsv(F("cog_categories"))
    if ch_ and cr_:
        fparts.append("<h3>" + L("COG 功能分类（eggNOG-mapper）", "COG functional categories (eggNOG-mapper)") + "</h3>" +
                      bar_figure([f"{x[0]}: {x[1]}" for x in cr_], [int(x[2]) for x in cr_], "genes", figdir,
                                 "cog_categories", (L("COG 功能分类。", "COG functional categories."),
                                                    L("eggNOG-mapper 给出的 COG 类别的基因数；一个基因可属于多个类别，S = 功能未知。",
                                                      "Genes per COG category from eggNOG-mapper; a gene can belong to "
                                                      "several categories, S = function unknown.")),
                                 wrap=48, values_text="{:.0f}") + table(ch_, cr_))
    lh, lr = read_tsv(F("eggnog_levels"))
    if lh and lr:
        fparts.append("<h3>" + L("eggNOG 注释所在的分类层级", "Taxonomic levels of the eggNOG annotations") + "</h3>" +
                      table(lh, lr))
    S.append(("function", L("06_Function 功能注释", "06_Function Functional annotation"), "".join(fparts)))
    # appendix
    vers = []
    for f in sorted(glob.glob(O("*/logs/versions/*.txt"))):
        txt = open(f, errors="replace").read().strip().splitlines()
        vers.append([os.path.relpath(f, a.outdir), (txt[0] if txt else "")[:120]])
    scripts = [[os.path.relpath(f, a.outdir)] for f in sorted(glob.glob(O("*/scripts/*.sh")))]
    cfg = html.escape(open(a.config).read()) if os.path.exists(a.config) else ""
    S.append(("appendix", L("附录：软件版本、脚本索引、参数快照", "Appendix: software versions, scripts, configuration"),
              "<h3>" + L("软件版本", "Software versions") + "</h3>" +
              table([L("文件", "File"), L("版本信息（首行）", "Version (first line)")], vers) +
              "<h3>" + L("各步骤脚本", "Scripts of all steps") + "</h3>" + table([L("脚本", "Script")], scripts) +
              f"<h3>config</h3><pre>{cfg}</pre>"))

    # the comparison covers EVM, refine and final: show it after the final result, before the functional annotation
    keys = [k for k, _, _ in S]
    if "qc" in keys and "finalset" in keys:
        qc = S.pop(keys.index("qc"))
        S.insert([k for k, _, _ in S].index("finalset") + 1, qc)

    res_html = resources_html(a.outdir, a.requests, figdir)
    if res_html:
        S.append(("resources", L("运行时间与资源", "Runtime and resources"), res_html))

    banner = ""
    if a.no_rna:   # sections that need RNA-seq are not produced; say so prominently
        S = [x for x in S if x[0] not in ("rnaqc", "rnaalign", "eviann")]
        banner = ("<div class='note'><b>" + L("无 RNA 模式", "No-RNA mode") + "</b>" +
                  L("（rna.allow_none: true）：未使用 RNA-seq。"
                    "BRAKER4 为 EP 模式（只用蛋白证据），不运行 EviAnn、StringTie，GeMoMa 不使用 RNA 剪接证据；"
                    "EVM 只整合 BRAKER4、GeMoMa 和 miniprot，权重调优只用蛋白证据打分；修正步骤没有 EviAnn 可比较，"
                    "所有位点按 R5 检查蛋白 / BRAKER4 证据；最终过滤以蛋白证据代替 RNA 证据。"
                    "结果没有 UTR、可变剪接和 lncRNA，基因结构的可靠性低于有 RNA 的注释。",
                    " (rna.allow_none: true): no RNA-seq was used. BRAKER4 runs in EP mode (protein evidence only), "
                    "EviAnn and StringTie are not run, GeMoMa uses no RNA splice evidence; EVM combines BRAKER4, GeMoMa "
                    "and miniprot only, and the weight tuning scores with protein evidence only; the refinement has no "
                    "EviAnn to compare with and checks every locus by R5 for protein / BRAKER4 evidence; the final "
                    "filter uses protein evidence instead of RNA evidence. The result has no UTRs, alternative isoforms "
                    "or lncRNAs, and gene structures are less reliable than with RNA.") + "</div>")

    def nav_item(k, t):
        step, _, rest = t.partition(" ")
        if re.match(r"^\d\d_", step) and rest:                  # "05_EVM 修正（…）" -> step label + title
            return f"<a href='#{k}'><span class='step'>{html.escape(step)}</span>{html.escape(rest)}</a>"
        return f"<a href='#{k}'>{html.escape(t)}</a>"
    other = a.other_name.get(LANG)
    switch = (f"<div class='lang'><a href='{html.escape(other)}'>" + L("English version", "中文版") + "</a></div>"
              if other else "")
    nav = (switch + "<div class='navtitle'>" + L("目录", "Contents") + "</div>" +
           "".join(nav_item(k, t) for k, t, _ in S) +
           "<a class='top' href='#top'>" + L("↑ 回到顶部", "↑ Back to top") + "</a>")
    body = number_figures("".join(f"<section id='{k}'><h2>{html.escape(t)}</h2>{c}</section>" for k, t, c in S))
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    title = L(f"{pfx} 注释报告", f"{pfx} annotation report")
    h1 = L(f"{pfx} 基因组注释统计报告", f"{pfx} genome annotation report")
    gen = L("由 Concordia 基因组注释流程（Snakemake）自动生成",
            "Generated automatically by the Concordia genome annotation workflow (Snakemake)")
    with open(out, "w") as fh:
        fh.write(f"<!doctype html><html lang='{LANG}'><head><meta charset='utf-8'>"
                 f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
                 f"<title>{title}</title><style>{CSS}</style></head>"
                 f"<body><div class='layout'><nav class='side'>{nav}</nav>"
                 f"<main><h1 id='top'>{h1}</h1>"
                 f"<p class='muted'>{gen}</p>{banner}{body}</main></div>"
                 f"<script>{NAV_JS}</script></body></html>")


def main():
    global LANG
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="pipeline result root (contains 00_Data ...)")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--braker4-mode", required=True)
    ap.add_argument("--eviann-tag", required=True, help="EviAnn input summary, e.g. '6 x bam + 1 x bam_isoseq'")
    ap.add_argument("--config", required=True)
    ap.add_argument("--eviann-proteins", default="homologs",
                    help="protein set given to EviAnn: swissprot | homologs | none (or free text)")
    ap.add_argument("--homologs", nargs="*", default=[], help="homologous species names (per-species QC dirs)")
    ap.add_argument("--no-rna", action="store_true", help="annotation without RNA-seq (rna.allow_none)")
    ap.add_argument("--requests", default="", help="07_Report/resource_requests.tsv (requested resources per rule)")
    ap.add_argument("--out-zh", help="Chinese report")
    ap.add_argument("--out-en", help="English report")
    a = ap.parse_args()
    outs = [(lang, p) for lang, p in (("zh", a.out_zh), ("en", a.out_en)) if p]
    if not outs:
        ap.error("give --out-zh and/or --out-en")
    # link from each version to the other (same directory: file name only)
    a.other_name = ({"zh": os.path.basename(a.out_en), "en": os.path.basename(a.out_zh)}
                    if a.out_zh and a.out_en else {})
    for lang, out in outs:
        LANG = lang
        build(a, out)


if __name__ == "__main__":
    main()
