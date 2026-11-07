#!/usr/bin/env python3
"""07_Report: build a self-contained HTML report (text + tables + figures) from all steps.

Skeleton (v0.3): genome / renaming, repeat summary, BRAKER4 mode, QC comparison (BUSCO figure),
refine decisions, filtering, function annotation rates, software versions and script index.
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
        return "<p class='muted'>（无数据）</p>"
    h = "".join(f"<th>{html.escape(c.lstrip('#'))}</th>" for c in header)
    b = "".join("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in r) + "</tr>" for r in rows[:limit])
    more = f"<p class='muted'>仅显示前 {limit} 行，共 {len(rows)} 行。</p>" if len(rows) > limit else ""
    return f"<div class='tbl'><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>{more}"


def fig_to_img(fig, figdir, name, caption=""):
    """Embed a figure; caption = (title, explanation) or text, shown as "图 N. title explanation"."""
    os.makedirs(figdir, exist_ok=True)
    fig.savefig(os.path.join(figdir, name + ".png"), dpi=150, bbox_inches="tight")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    img = f"<img alt='{name}' src='data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}'>"
    title, text = caption if isinstance(caption, tuple) else ("", caption)
    cap = (f"<figcaption><b>图 {FIG_MARK}. {html.escape(title)}</b> {html.escape(text)} "
           f"<span class='muted'>（{html.escape(name)}.png）</span></figcaption>")
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
        "BUSCO 评估。", "每行一个注释来源（或基因组）× 一个 BUSCO 数据集；条形按单拷贝完整（S）、重复完整（D）、"
        "片段（F）、缺失（M）的百分比堆叠，右侧为完整度 C = S + D。"))


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
    iso_txt = ("左：每个基因的转录本数（1、2、3、4、≥5 个）占基因的百分比；" if iso_sets else
               "（每个基因只有一个转录本，不画转录本数分布。）")
    return fig_to_img(f, figdir, name, (
        f"{what}基因结构分布。",
        iso_txt + "中：每个转录本的外显子数；右：CDS 长度（氨基酸）。后两幅为各注释集内转录本的百分比，"
        "便于比较规模不同的注释集；超过上限的值计入最后一格。"))


REF_COMPARE_HELP = """
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

TUNE_HELP = """
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
    rules = "".join([
        row("R1", "1 个 EVM 基因 : 1 个 EviAnn 基因，EviAnn 有一条转录本的 CDS 内含子链与 EVM 完全相同",
            "保留 EVM 的 CDS；外显子和 UTR 取自这条 EviAnn 转录本；EviAnn 的其它完整转录本（完整 ORF ≥ min_cds_aa）作为异构体加入"
            f"（import_utr_isoforms = {v['import_utr_isoforms']}）。若 EVM 的 CDS 不落在 EviAnn 外显子内，只保留 EVM。不依赖 RNA 支持"),
        row("R2", "只有 EviAnn 基因，没有重叠的 EVM 基因",
            f"同时满足才加入：① 有完整 ORF ≥ {v['min_cds_aa']} aa；② 有 RNA 证据（与 StringTie 转录本重叠，或"
            f"按下表“被支持”）；③ CDS 与 TE 重叠 ≤ {v['max_te_overlap']}。否则跳过（skip）"),
        row("R3", "1 : 1，但 CDS 内含子链不同",
            "EviAnn 有被支持的完整转录本、EVM 不被支持 → 用 EviAnn 替换（replace_with_eviann）；两者都被支持 → 保留 EVM，"
            "并把被支持的 EviAnn 转录本加为异构体（keep_evm+eviann_isoforms）；其它 → 保留 EVM"),
        row("R4 拆分", "1 个 EVM 基因 : 多个 EviAnn 基因",
            "每个 EviAnn 基因都有被支持的转录本、且没有 reads 把它们连起来 → 用这些 EviAnn 基因替换；有 reads 连接 → 保留 EVM"
            "（说明确实是一个基因）；其它 → 保留 EVM，标记 refine_support=needs_review"),
        row("R4 融合", "多个 EVM 基因 : 1 个 EviAnn 基因",
            "EviAnn 有被支持、且 CDS 覆盖所有这些 EVM 基因的转录本 → 用 EviAnn 替换；否则保留 EVM，标记 needs_review"),
        row("R4 复杂", "多 : 多", "保留 EVM，标记 needs_review"),
        row("R5", "只有 EVM 基因，没有重叠的 EviAnn 编码基因",
            "有以下任一证据则保留并记录在 refine_support：蛋白（miniprot / GeMoMa 的 CDS 同链重叠）、RNA（StringTie 重叠或“被支持”）、"
            "BRAKER4 gene_support 的 RNA 支持；都没有 → "
            + ("删除" if v["remove_unsupported_evm"] == "true" else "保留并标记 refine_support=none")
            + f"（remove_unsupported_evm = {v['remove_unsupported_evm']}）"),
        row("R6", "EviAnn 的 lncRNA", "在更早的步骤已分出，单独输出为 lncRNA 文件，不参与修正"),
    ])
    support = "".join([
        row("有长读", "多外显子模型",
            "长读在模型范围内的内含子与模型的内含子完全一致；Iso-Seq HQ/LQ 转录本按 full_length_coverage（全长 reads 数）加权，"
            f"合计 ≥ min_lr_support = {v['min_lr_support']}"),
        row("有长读", "单外显子模型",
            f"无剪接的长读覆盖模型长度 ≥ single_exon_overlap = {v['single_exon_overlap']}（加权计数同上）"),
        row("只有短读（B）", "1 个内含子",
            f"该内含子的剪接 reads ≥ min_sr_junction_reads = {v['min_sr_junction_reads']}，且 ≥ sr_min_junction_ratio = "
            f"{v['sr_min_junction_ratio']} × 与它重叠的最强竞争内含子的 reads 数（排除低频的错误剪接）"),
        row("只有短读（A 且 B）", "≥ 2 个内含子",
            "每个内含子都满足 B；并且有一条 StringTie 转录本（短读组装，cov ≥ sr_stringtie_min_cov = "
            f"{v['sr_stringtie_min_cov']}，同链或无链）在模型范围内的内含子链与模型完全一致，证明这些内含子属于同一条转录本"),
        row("只有短读（C）", "单外显子模型",
            f"有单外显子 StringTie 转录本（cov 达标）覆盖模型长度 ≥ {v['single_exon_overlap']}，且模型范围内没有达到 B 绝对门槛的剪接位点"),
        row("无 RNA 数据", "—", "一律视为不被支持"),
    ])
    th = lambda *c: "<thead><tr>" + "".join(f"<th>{x}</th>" for x in c) + "</tr></thead>"
    return ("<details><summary><b>修正规则详解（点击展开）</b></summary>"
            "<p><b>位点与比较方式：</b>EVM 基因和 EviAnn 编码基因若在同一条链上 CDS 有重叠，归入同一个位点（可连锁合并），"
            "规则按位点内两边的基因数选择。模型之间比较 CDS 的内含子链（EVM 模型没有 UTR，比较完整外显子很难一致）。"
            f"EviAnn 转录本须有完整 ORF（起始和终止密码子、无内部终止）且 ≥ min_cds_aa = {v['min_cds_aa']} aa 才会被采用。</p>"
            "<h4>规则</h4><table>" + th("规则", "情形", "处理") + "<tbody>" + rules + "</tbody></table>"
            "<h4>“被支持”的判定（RNA 支持）</h4><table>" + th("数据", "模型", "条件") + "<tbody>" + support + "</tbody></table>"
            "<p><b>跨基因连接（R4）：</b>有长读时为同时含两个相邻基因外显子的长读；只有短读时为一端在基因 A 外显子、另一端在基因 B 外显子的剪接 reads。</p>"
            "<p><b>输出记录：</b>基因属性 refine_rule（适用的规则）、eviann_support（是否与 EviAnn 一致或来自 EviAnn）、refine_support"
            "（R5 的证据或 needs_review）；05_EVM/refine/refine_decisions.tsv 每个位点一行，support_basis 列记录支持依据（long / A+B / B / C）。"
            "表中数值为本次运行的配置（refine 段）。</p></details>")


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
        note = ("<div class='note'><b>注意：</b>Iso-Seq 样本 " + html.escape("、".join(warn)) + " 没有任何一条记录带全长 reads 数"
                "（表头没有 full_length_coverage，样本表也没有 counts 列），每条记录按 1 计数。若这些记录是原始 reads（CCS / FLNC），"
                "这是正确的；若是聚类后的 HQ 转录本，长读支持会被明显低估（修正 R3/R4 偏保守、调优 P 偏低）。"
                "请在 rna_long.tsv 的 counts 列填写 isoseq cluster 的 cluster_report.csv 或 isoseq collapse 的 flnc_count.txt，"
                "或将 refine.min_lr_support 设为 1。</div>")
    return ("<h3>长读的全长 reads 权重</h3><p class='muted'>修正步骤（R3/R4）和权重调优按\"全长 reads\"计算长读支持："
            "每条比对上的记录的权重取自样本表 counts 列的计数文件，或 FASTA/FASTQ 表头的 full_length_coverage，都没有时按 1 计。</p>" +
            note + table(h, r))


def tune_html(path, best, figdir):
    """EVM weight-tuning results: score table, F bar chart, the weights finally used."""
    if not os.path.exists(path):
        return "<p class='muted'>未进行权重调优（evm.tune.enabled: false），使用配置中的权重。</p>"
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
            "EVM 权重调优。", "每行一组权重（c00 为配置中的基线），按 F 从高到低排列；P = 多外显子模型的内含子链被任一证据完全确认的比例，"
            "R = 证据位点被 EVM 模型覆盖的比例，F = 2PR/(P+R)。横轴只显示数值所在区间，便于看出差别。"))
    used = [l.rstrip("\n").split("\t") for l in open(best) if l.strip()] if os.path.exists(best) else []
    base = next((x for x in r if x[1] == "c00"), None)
    cmp_ = (f"最优组合 <b>{html.escape(r[0][1])}</b>（F = {r[0][h.index('F')]}）；配置中的基线权重 c00 排第 "
            f"{base[0]} 位（F = {base[h.index('F')]}）。" if base else "")
    return ("<h3>权重调优</h3>" + TUNE_HELP + "<p class='muted'>" + html.escape("; ".join(note)) + "</p><div class='note'>" + cmp_ +
            "</div>" + fig + table(h, r, limit=40) + "<p>全基因组 EVM 实际使用的权重：</p>" +
            table(["类别", "来源", "权重"], used))


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
        rows.append([code, ", ".join(seqs[:3]) + (f" 等 {len(seqs)} 条" if len(seqs) > 3 else ""), v["G"], v["L"],
                     first.get((code, "G"), "-")])
    return ("<h3>基因命名与各染色体基因数</h3><p class='muted'>编码基因 &lt;前缀&gt;&lt;染色体号&gt;G&lt;编号&gt;，"
            "lncRNA 用 L；未挂载到染色体的序列代码为 00；新旧 ID 对照见 05_EVM/final/id_map.tsv。</p>" +
            table(["染色体代码", "序列", "编码基因", "lncRNA", "第一个基因"], rows))


def source_support_html(support_tsv, combos_tsv, figdir, name):
    sh_, sr_ = read_tsv(support_tsv)
    ch_, cr_ = read_tsv(combos_tsv)
    if not sh_:
        return ""
    fig = ""
    if cr_:
        top = cr_[:12]
        fig = bar_figure([x[0] for x in top], [int(x[1]) for x in top], "transcripts", figdir,
                         f"{name}_source_combinations", (
                             ("EVM 结果" if name == "evm" else "最终结果") + "的来源组合。",
                             "每个转录本按哪些输入来源（BRAKER4 / GeMoMa / EviAnn）有 CDS 结构完全相同的模型分组；"
                             "no exact match = 与任何来源都不完全一致（EVM 拼合或修正产生的结构）。"), values_text="{:.0f}")
    return ("<h3>与各输入来源的一致性</h3><p class='muted'>exact = 该来源有 CDS 结构完全相同的转录本；"
            "overlap = 与该来源的 CDS 在同一条链上有重叠。组合 = 有完全一致模型的来源组合。</p>" +
            table(sh_, sr_) + fig + table(ch_, cr_, limit=30))


SECTION_WHAT = {"eviann": "EviAnn 注释（编码基因与 lncRNA）", "braker4": "BRAKER4 注释", "evm": "EVM 结果",
                "gemoma": "同源注释（各物种单独结果与合并结果）", "final": "最终结果（主结果）"}


def annotation_section(O, name, stats_tsv, dist_tsv, qc_dir, figdir, intro, extra=""):
    """Statistics table + distributions + protein-set quality + AGAT for one annotation source."""
    h, r = read_tsv(stats_tsv)
    _, dist = read_tsv(dist_tsv)
    parts = [intro, "<h3>基因结构统计</h3>",
             "<p class='muted'>UTR 优先取 GFF 中的 UTR 特征，否则按外显子超出 CDS 的部分判断（区分正负链）。</p>",
             table(h, r), distribution_figure(dist, figdir, f"{name}_distributions", SECTION_WHAT.get(name, ""))]
    q = protein_qc_rows(qc_dir)
    if q:
        parts.append("<h3>蛋白集质量（最长异构体）</h3>" + table(["评估", "结果"], q))
    parts.append(extra)
    agat_f = os.path.join(qc_dir, "agat", "stats.txt")
    if os.path.exists(agat_f):
        parts.append("<details><summary>AGAT 统计全文</summary><pre>" +
                     html.escape(open(agat_f).read()[:20000]) + "</pre></details>")
    return "".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="pipeline result root (contains 00_Data ...)")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--braker4-mode", required=True)
    ap.add_argument("--eviann-tag", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--eviann-proteins", default="同源物种蛋白", help="protein set given to EviAnn")
    ap.add_argument("--homologs", nargs="*", default=[], help="homologous species names (per-species QC dirs)")
    ap.add_argument("--no-rna", action="store_true", help="annotation without RNA-seq (rna.allow_none)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    O = lambda *p: os.path.join(a.outdir, *p)
    figdir = os.path.join(os.path.dirname(a.out), "figures")
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
        hom_html = ("<h3>同源物种注释预处理</h3><p class='muted'>只保留蛋白编码基因（及其有 CDS 的转录本），"
                    "去掉转座子基因、假基因、ncRNA 基因等；GeMoMa 和蛋白库都使用预处理后的注释。</p>" +
                    table(["项目"] + list(per), [[k] + [per[n].get(k, "0") for n in per] for k in items]))
    S.append(("data", "00_Data 输入与序列名检查", table(["项目", "值"], rn) + hom_html))

    # 00_Data: genome assembly quality (every assessment result appears here)
    AQ = lambda *p: O("00_Data/qc/assembly", *p)
    qh, qr = read_tsv(O("00_Data/qc/assembly_quality.tsv"))
    parts = ["<p>对用户提供的原始组装（原始序列名）进行评估。连续性：scaffold / contig 两级，contig 在连续 N 处拆分；"
             "完整度：BUSCO 基因组模式、Merqury k-mer 完整度（可选）、LAI（EDTA 完成后计算）；准确度：Merqury QV（可选）；"
             "端粒：tidk；污染：FCS-GX（可选，只报告不删除）。</p>",
             "<h3>汇总</h3>" + table(qh, qr)]
    skipped = [x[1] for x in qr if len(x) > 2 and x[2] == "not run"]
    if skipped:
        parts.append("<div class='note'>未运行（未配置所需数据或数据库）：" + html.escape("、".join(skipped)) + "</div>")

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
            "序列累计长度。", "序列按长度从长到短排列（横轴为对数刻度），纵轴为累计长度占组装总长的百分比；"
            "虚线 50% 和点线 90% 与曲线的交点分别对应 L50 和 L90。"))
    sk_h, sk_r = read_tsv(O("00_Data/qc/seqkit/genome.stats.tsv"))
    parts.append("<h3>连续性</h3>" + table(sh_, sr_) + fig +
                 "<p class='muted'>seqkit stats -a（重命名和长度过滤后的内部基因组）：</p>" + table(sk_h, sk_r))
    if seqs:
        parts.append("<h3>各序列</h3>" + table(["seq_id", "length", "GC_percent", "N_bases", "gaps", "gap_length",
                                                "contigs"], seqs, limit=60))

    # BUSCO genome mode
    brows = []
    for f_ in sorted(glob.glob(O("00_Data/qc/busco_genome/*/short_summary.txt"))):
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", open(f_).read())
        if m:
            brows.append(["genome", os.path.basename(os.path.dirname(f_)), *m.groups()])
    if brows:
        bh = ["source", "lineage", "C", "S", "D", "F", "M", "n"]
        parts.append("<h3>BUSCO（基因组模式）</h3>" + busco_figure(bh, brows, figdir, "busco_genome", (
            "基因组 BUSCO（genome 模式）。", "每行一个 BUSCO 数据集；按单拷贝完整（S）、重复完整（D）、片段（F）、缺失（M）的百分比堆叠，"
            "右侧为完整度 C = S + D。这是基因组组装完整性的基线，注释的 BUSCO 一般不会高于它。")) + table(bh, brows))

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
        lai_help = ("<p class='muted'>LAI（LTR Assembly Index，Ou et al. 2018, NAR 46:e126）用完整 LTR 反转录转座子占全部 LTR 序列的比例"
                    "衡量组装对重复区的连续性：raw_LAI = 完整 LTR 长度 / 全部 LTR 长度 × 100；LAI 再按 LTR 序列一致性校正，"
                    "排除近期大量扩增的家族造成的偏高。参考分级：draft &lt; 10，reference 10–20，gold ≥ 20。"
                    "完整 LTR 来自 EDTA，全部 LTR 来自 RepeatMasker（EDTA 库）。</p>")
        first = open(lai_f).readline()
        if first.startswith("# LAI not applicable"):
            parts.append("<h3>LAI</h3>" + lai_help + "<div class='note'>LAI 不适用：" +
                         html.escape(first.split(":", 1)[1].strip().rstrip(".")) +
                         "。LTR 反转录转座子占基因组不足 5% 时，LAI 无法可靠估计（单条染色体测试或 LTR 很少的物种常见）。</div>")
        else:
            lh, lr = read_tsv_ws(lai_f)
            whole = [x for x in lr if x and x[0] == "whole_genome"]
            parts.append("<h3>LAI</h3>" + lai_help +
                         (f"<div class='note'>全基因组 LAI = <b>{html.escape(whole[0][-1])}</b>"
                          f"（raw_LAI = {html.escape(whole[0][-2])}）</div>" if whole else "") +
                         "<details><summary>按序列 / 窗口的 LAI</summary>" + table(lh, lr, limit=200) + "</details>")

    # Telomeres
    ts = kv_file(AQ("telomere", "telomere_summary.tsv"))
    th, tr = read_tsv(AQ("telomere", "telomere_per_sequence.tsv"))
    if ts or th:
        parts.append("<h3>端粒（tidk）</h3>" + table(["项目", "值"], ts) + table(th, tr, limit=60))

    # FCS-GX
    fh_, fr_ = read_tsv(AQ("fcs_gx", "fcs_gx_report.txt"))
    if fh_:
        parts.append("<h3>污染筛查（FCS-GX）</h3><p>只报告，流程不会自动删除任何序列。</p>" +
                     (table(fh_, fr_) if fr_ else "<p>未发现污染。</p>"))

    S.append(("assembly", "00_Data 基因组组装质量", "".join(parts)))

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
    S.append(("rnaqc", "00_Data RNA-seq 原始数据",
              "<h3>短读（fastp）</h3>" + table(["sample", "reads_before", "reads_after", "Q30_%", "GC_%"], fastp) +
              "<h3>长读（NanoPlot）</h3>" + table(["sample", "reads", "mean_length", "N50", "mean_quality"], nano)))
    # 01_Repeat
    R = lambda n: O("01_Repeat/stats", n + ".tsv")
    rparts = ["<p>EDTA 构建物种特异 TE 库并做全基因组 TE 注释；RepeatMasker 用该库对基因组软屏蔽"
              "（同时屏蔽简单重复和低复杂度序列）。重叠的注释只计一次。</p>"]
    h, r = read_tsv(R("overview"))
    if h:
        rparts.append("<h3>总览</h3>" + table(h, r))
    h, r = read_tsv(R("te_classes"))
    if h and r:
        tops = [x for x in r if x[1] == "(all)"]
        fig = bar_figure([x[0] for x in tops], [float(x[4]) if x[4] != "-" else 0 for x in tops], "% of genome",
                         figdir, "te_classes", ("各类转座子占基因组的比例。",
                                                "EDTA 注释中各 TE 大类（各超家族合计）覆盖的碱基占基因组长度的百分比；重叠部分只计一次。"),
                         values_text="{:.2f}%")
        rparts.append("<h3>TE 类别组成</h3>" + fig + table(h, r))
    h, r = read_tsv(R("te_methods"))
    if h and r:
        rparts.append("<h3>结构鉴定（完整元件）与同源注释</h3><p class='muted'>structural = 由结构特征鉴定的完整元件；"
                      "homology = 根据 TE 库同源比对注释的片段。</p>" + table(h, r))
    h, r = read_tsv(R("library"))
    if h and r:
        rparts.append("<h3>TE 库</h3><p class='muted'>共 " + str(sum(int(x[1]) for x in r)) + " 条一致序列。</p>" +
                      table(h, r))
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
        ax.legend(fontsize=8, ncol=4, frameon=False)
        rparts.append("<h3>重复序列景观（Kimura 距离）</h3><p class='muted'>距离越小，插入越年轻。</p>" +
                      fig_to_img(f, figdir, "repeat_landscape", (
                          "重复序列景观。", "横轴为拷贝与其一致序列之间的 Kimura 距离（经 CpG 校正，%），纵轴为该距离区间内各 TE 类别占基因组的百分比。"
                          "距离越小代表插入越年轻；0 附近的峰常来自近期活跃或刚扩增的家族。")))
    age_f = R("ltr_age")
    if os.path.exists(age_f):
        lines = [l.rstrip("\n").split("\t") for l in open(age_f) if l.strip()]
        note = [l[0].lstrip("# ") for l in lines if l[0].startswith("#")]
        body = [l for l in lines[1:] if not l[0].startswith("#")][1:] if len(lines) > 2 else []
        if body:
            by_sf = {}
            for x in body:
                by_sf.setdefault(x[1], []).append(float(x[3]))
            summ = [[sf, len(v), f"{sorted(v)[len(v) // 2]:.2f}", f"{min(v):.2f}", f"{max(v):.2f}"]
                    for sf, v in sorted(by_sf.items(), key=lambda t: -len(t[1]))]
            fig = ""
            n_ltr = sum(len(v) for v in by_sf.values())
            if plt is not None and n_ltr < 20:
                fig = f"<p class='muted'>只有 {n_ltr} 个完整 LTR，数量太少，不画分布图，见下表。</p>"
            elif plt is not None:
                f, ax = plt.subplots(figsize=(7, 3.2))
                for i, (sf, v) in enumerate(sorted(by_sf.items(), key=lambda t: -len(t[1]))[:4]):
                    ax.hist(v, bins=40, alpha=0.6, label=sf, color=["#2f6fb0", "#e0a43a", "#5a9e6f", "#c8553d"][i])
                ax.set_xlabel("insertion time (Mya)")
                ax.set_ylabel("intact LTRs")
                ax.legend(frameon=False)
                fig = fig_to_img(f, figdir, "ltr_insertion_time", (
                    "完整 LTR 反转录转座子的插入时间。", "按超家族（最多 4 个）分别统计；插入时间 T = K / (2u)，K 为两端 LTR 之间的分化度，"
                    "u 为替换速率（EDTA -u，见配置）。"))
            rparts.append("<h3>完整 LTR 插入时间</h3><p class='muted'>" + html.escape("; ".join(note)) + "</p>" + fig +
                          table(["超家族", "完整 LTR 数", "中位数 (Mya)", "最小", "最大"], summ))
        else:
            rparts.append("<h3>完整 LTR 插入时间</h3><p class='muted'>未找到完整 LTR 列表或其中没有插入时间。</p>")
    sumf = glob.glob(O("01_Repeat/EDTA/*.TEanno.sum"))
    if sumf:
        rparts.append("<details><summary>EDTA 注释汇总原文</summary><pre>" +
                      html.escape(open(sumf[0]).read()[:20000]) + "</pre></details>")
    S.append(("repeat", "01_Repeat 重复序列", "".join(rparts)))
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
            "RNA-seq 比对率。", "每行一个样本（括号内为平台）。蓝色 = 比对上的 reads 占总 reads 的百分比；"
            "橙色 = 比对上的 reads 中含剪接（CIGAR 中有 N）的比例，长读通常远高于短读。"))
    S.append(("rnaalign", "02_RNA RNA-seq 比对",
              "<p>短读用 HISAT2 比对，长读用 minimap2 比对（预设按平台选择）。比对率 = 比对上的 reads / 总 reads"
              "（samtools stats，不含 secondary/supplementary）；unique/multi 来自 HISAT2 汇总；"
              "spliced % = 含内含子（CIGAR 中有 N）的主比对占比对上 reads 的比例；error rate 为错配碱基 / 比对碱基。</p>"
              "<h3>各样本</h3>" + fig + table(ah, ar) +
              "<h3>合并 BAM 及其下游用途</h3>" + table(mh, mr) + long_read_weights_html(O("02_RNA/qc/long_read_weights.tsv"))))
    # 02-04: annotation sources (same layout for EviAnn, BRAKER4, GeMoMa)
    S.append(("eviann", "02_RNA EviAnn 注释", annotation_section(
        O, "eviann", O("02_RNA/EviAnn/stats/eviann_stats.tsv"), O("02_RNA/EviAnn/stats/eviann_distributions.tsv"),
        O("02_RNA/EviAnn/qc"), figdir,
        f"<p>EviAnn 以 RNA 比对（reads 标签 <b>{html.escape(a.eviann_tag)}</b>）为证据进行注释；确定 CDS 所用的蛋白集："
        f"<b>{html.escape(a.eviann_proteins)}</b>。"
        "结果按是否含 CDS 拆分为编码基因（进入 EVM）和 lncRNA（单独输出）。</p>")))

    sp_h, sp_r = read_tsv(O("03_Ab_initio/stats/braker4_support.tsv"))
    sp_note = ""
    sp_file = O("03_Ab_initio/stats/braker4_support.tsv")
    if os.path.exists(sp_file):
        sp_note = "".join(f"<p class='muted'>{html.escape(l[1:].strip())}</p>" for l in open(sp_file) if l.startswith("#"))
    sp_r = [r for r in sp_r if not r[0].startswith("#")]
    S.append(("braker4", "03_Ab_initio BRAKER4 注释", annotation_section(
        O, "braker4", O("03_Ab_initio/stats/braker4_stats.tsv"), O("03_Ab_initio/stats/braker4_distributions.tsv"),
        O("03_Ab_initio/qc"), figdir,
        f"<div class='note'>BRAKER4 模式：<b>{html.escape(a.braker4_mode)}</b>（只使用 PacBio 长读；"
        "ONT 不进入 BRAKER4）。输入为软屏蔽基因组、OrthoDB 与同源物种蛋白。</div>",
        extra="<h3>证据支持（gene_support.tsv）</h3>" + sp_note +
              "<p class='muted'>BRAKER4 未说明该文件的列含义，此处逐列统计非零/为真的模型比例。</p>" +
              table(sp_h, sp_r))))

    rf_h, rf_r = read_tsv(O("04_Homology/stats/gemoma_references.tsv"))
    at_h, at_r = read_tsv(O("04_Homology/stats/gemoma_attributes.tsv"))
    mp_h, mp_r = read_tsv(O("04_Homology/stats/miniprot_species.tsv"))
    sp_qc = []
    for h_ in a.homologs:
        for k, v in protein_qc_rows(O("04_Homology/qc", f"GeMoMa_{h_}")):
            sp_qc.append([h_, k, v])
    S.append(("gemoma", "04_Homology 同源注释（GeMoMa / miniprot）", annotation_section(
        O, "gemoma", O("04_Homology/stats/gemoma_stats.tsv"), O("04_Homology/stats/gemoma_distributions.tsv"),
        O("04_Homology/qc/GeMoMa_merged"), figdir,
        "<p>每个同源物种单独运行一次 GeMoMa：把该物种的基因模型投射到目标基因组，并用 RNA-seq 的剪接位点和覆盖度"
        "（ERE，只提取一次）校正。各物种的预测再用 GeMoMa GAF + AnnotationFinalizer 合并为 GeMoMa_merged，"
        "作为 EVM 的输入。下表中 GeMoMa_&lt;物种&gt; 为各物种单独的结果；蛋白集质量一栏为合并结果。</p>",
        extra="<h3>各物种单独结果的蛋白集质量</h3>" + (table(["物种", "评估", "结果"], sp_qc) if sp_qc else
                                                      "<p class='muted'>未做按物种的评估（qc.per_homolog: false）。</p>") +
              "<h3>各物种对合并结果的贡献</h3><p class='muted'>merged_representative = 合并结果中代表预测来自该物种的转录本；"
              "merged_supported = 该物种也预测出相同 CDS 的转录本（含 GAF alternative）。</p>" + table(rf_h, rf_r) +
              "<h3>合并结果的证据属性</h3><p class='muted'>evidence = 预测出该转录本的物种数；tie = 有 RNA 支持的内含子比例；"
              "tpc = 被 RNA 覆盖的转录本比例（后两项只在使用 RNA 证据时有意义）。</p>" + table(at_h, at_r) +
              "<h3>miniprot 蛋白比对（按物种）</h3><p class='muted'>每个物种的蛋白单独比对到基因组；合并后作为 EVM 的 "
              "PROTEIN 证据，也用于修正（R 规则）和权重调优。</p>" + table(mp_h, mp_r))))

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
                    "与参考注释的一致性（召回）。", "gene_recall = 参考基因中与该注释集有同链 CDS 重叠的比例；intron_recall = 参考的 CDS 内含子"
                    "被完全相同地预测出的比例；tx_recall = 参考转录本的完整 CDS 结构被完全相同地预测出的比例。精确度见下表。"
                    "参考注释本身也是预测结果时，反映的是一致程度而非准确率。"))
            ref_html = ("<h3>与参考注释的一致性（CDS 水平）</h3>" + REF_COMPARE_HELP + "<p class='muted'>" + html.escape("; ".join(note)) +
                        "</p>" + fig + table(rh, rr))
    ch, cr = read_tsv(O("05_EVM/stats/all_sets_stats.tsv"))
    S.append(("qc", "各注释来源对比（BUSCO / OMArk / 基因结构）",
              "<p>所有来源使用同一组用户指定的 BUSCO lineage。</p>" + busco_figure(h, r, figdir, "busco_comparison", (
                  "各注释来源的蛋白集 BUSCO。", "每行一个注释来源 × 一个 BUSCO 数据集，输入为各来源每个基因的最长蛋白；"
                  "按 S / D / F / M 堆叠，右侧为完整度 C。各来源使用同一组数据集，可以直接比较。")) + table(h, r) +
              "<h3>基因结构对比</h3><p class='muted'>EviAnn 一列只含编码基因；final 为过滤和重新命名后的最终基因集。</p>" +
              table(ch, cr) + ref_html))
    # 05_EVM: EVM integration
    impl = kv_file(O("05_EVM/logs/evm_impl_check.txt"))
    impl_txt = ""
    if impl:
        bits = [impl[0][0]] + impl[0][1].split("\t")
        impl_txt = (f"<div class='note'>EVM 实现：<b>{html.escape(bits[0])}</b>（{html.escape(bits[-1])}）"
                    f"<br><span class='muted'>{html.escape(bits[1] if len(bits) > 2 else '')}</span></div>")
    w_rows = [l.rstrip("\n").split("\t") for l in open(O("05_EVM/inputs/weights.txt"))
              if l.strip()] if os.path.exists(O("05_EVM/inputs/weights.txt")) else []
    S.append(("evm", "05_EVM EVM 整合", annotation_section(
        O, "evm", O("05_EVM/stats/evm_stats.tsv"), O("05_EVM/stats/evm_distributions.tsv"), O("05_EVM/qc/evm"),
        figdir, impl_txt + tune_html(O("05_EVM/tune/tuning_results.tsv"), O("05_EVM/tune/best_weights.txt"), figdir)
        + "<h3>输入权重</h3>" + table(["类别", "来源", "权重"], w_rows),
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
    rparts = ["<p>EVM 与 EviAnn 编码基因按 CDS 在同一条链上的重叠划分为位点，按 CDS 内含子链比较：</p><ul>"
              "<li>R1 一致：保留 EVM 的 CDS，引入 EviAnn 的 UTR 和其它完整异构体</li>"
              "<li>R2 EVM 缺失：EviAnn 基因有完整 ORF、RNA 支持且 CDS 与 TE 重叠不超过阈值时补入</li>"
              "<li>R3 结构冲突：EviAnn 内含子链被 RNA 支持而 EVM 不被支持时替换；两者都被支持时把 EviAnn 模型加为异构体</li>"
              "<li>R4 拆分 / 融合：有 RNA 支持时采用 EviAnn 的结构，有 reads 连接两个基因时保留融合；证据不足时保留 EVM 并标记 needs_review</li>"
              "<li>RNA 支持：有长读时为长读内含子链与模型完全一致；只有短读时，1 个内含子看剪接 reads（B：绝对数与相对竞争内含子的比例），"
              "≥ 2 个内含子还要求 StringTie 组装出相同的内含子链（A 且 B），单外显子模型看单外显子 StringTie 转录本（C）；"
              "决策表 support_basis 列记录依据</li>"
              "<li>R5 EVM 独有：检查蛋白、RNA、BRAKER4 证据，没有证据的标记 refine_support=none</li>"
              "<li>R6 lncRNA：单独输出，不参与合并</li></ul>" + refine_rules_html(a.config)]
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
            ("修正规则的处理结果。", "每个位点（EVM 与 EviAnn 编码基因按同链 CDS 重叠分组）按适用的规则和处理方式计数；"
             "merge_utr (+isoforms) 合并了引入不同数目异构体的情况，明细见下表。颜色区分规则。"),
            colors=cols, legend=[(k, v) for k, v in (("R1", "#2f6fb0"), ("R2", "#5a9e6f"), ("R3", "#e0a43a"),
                                                       ("R4", "#c8553d"), ("R5", "#8c6bb1"))
                                 if any(x[0].startswith(k) for x in gkeys)], values_text="{:.0f}"))
        rparts.append(table(["规则", "处理", "位点数"], [[k[0], k[1], c[k]] for k in keys]))
        review = [x for x in r if "needs_review" in x[h.index("note")]]
        if review:
            rparts.append(f"<h3>需要人工检查的位点（{len(review)} 个）</h3>" + table(h, review, limit=50))
        rparts.append("<details><summary>全部决策（前 200 行）</summary>" + table(h, r, limit=200) + "</details>")
    S.append(("refine", "05_EVM 修正（EVM vs EviAnn）", "".join(rparts)))

    # filtering
    h, r = read_tsv(O("05_EVM/final/removed_genes.tsv"))
    fparts = []
    if h:
        c = Counter(x[1].split("_")[0] for x in r)
        fparts.append(table(["原因", "基因数"], sorted(c.items())))
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
            summ = [["参与过滤的基因", len(body)],
                    [f"TE 重叠 > {thr}", len(high_te)],
                    ["其中因有 RNA 证据而保留", len(rescued)],
                    ["　StringTie 支持", by["rna_stringtie"]], ["　EviAnn 支持", by["rna_eviann"]],
                    ["　BRAKER4 gene_support 支持", by["rna_braker4"]],
                    ["有 TEsorter 命中（仅记录）", sum(x[ix["tesorter"]] != "-" for x in body)]]
            fparts.append("<h3>TE 过滤明细</h3>" + "".join(f"<p class='muted'>{html.escape(n_)}</p>" for n_ in notes) +
                          table(["项目", "数目"], summ))
    S.append(("filter", "05_EVM 过滤",
              "<p>基因 CDS 与 TE 重叠超过阈值且没有 RNA 证据时才剔除；CDS 过短的模型单独剔除；TEsorter 只记录。</p>" +
              "".join(fparts)))

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
            kind = "序列"
        else:
            n_rec = sum(1 for l in open(f_) if "\tgene\t" in l)
            kind = "基因"
        files.append([os.path.relpath(f_, a.outdir), f"{size / 1e6:.2f} MB", f"{n_rec} {kind}"])
    S.append(("finalset", "05_EVM 最终注释结果", annotation_section(
        O, "final", O("05_EVM/stats/final_stats.tsv"), O("05_EVM/stats/final_distributions.tsv"),
        O("05_EVM/qc/final"), figdir,
        "<p>最终基因集：EVM 整合 → 与 EviAnn 比较修正 → TE / 长度过滤 → 重新命名并还原原始序列名。"
        "lncRNA 来自 EviAnn，单独输出。</p>"
        "<div class='note'>主结果（<code>" + html.escape(a.prefix) + ".gff3</code> 及 cds/pep/cdna）每个基因只保留一个转录本："
        "CDS 最长者，CDS 等长时取转录本最长者；本节统计、BUSCO/OMArk 和功能注释都基于主结果。"
        "包含全部异构体的结果另存为 <code>" + html.escape(a.prefix) + ".all_isoforms.*</code>。</div>" +
        (table(["项目", "数目"], [[k, v] for k, v in kv_file(O("05_EVM/final/isoform_selection.tsv"))
                                 if k != "item"] +
               ([["transcripts_with_internal_stop (kept, internal_stop=yes)",
                  str(max(0, len(open(O("05_EVM/final/internal_stops.tsv")).read().splitlines()) - 1))]]
                if os.path.exists(O("05_EVM/final/internal_stops.tsv")) else []))
         if os.path.exists(O("05_EVM/final/isoform_selection.tsv")) else "") +
        "<h3>各阶段基因数</h3>" + table(["阶段", "基因", "转录本"], prog),
        extra=source_support_html(O("05_EVM/stats/final_source_support.tsv"),
                                  O("05_EVM/stats/final_source_combinations.tsv"), figdir, "final") +
              naming_html(O("05_EVM/final/id_map.tsv"), a.prefix) +
              "<h3>输出文件</h3>" + table(["文件", "大小", "记录数"], files))))

    # function
    F = lambda n: O("06_Function/stats", n + ".tsv")
    fparts = ["<p>对最终基因集主结果的蛋白（每个基因一个转录本）进行 InterProScan 和 eggNOG-mapper 注释，"
              "统计均以基因为单位。</p>"]
    oh, orows = read_tsv(F("overview"))
    if oh:
        keep = [x for x in orows if x[0] in ("InterProScan only", "both", "eggNOG-mapper only", "neither (unannotated)")]
        fig = bar_figure([x[0] for x in keep], [float(x[2]) for x in keep], "% genes", figdir, "function_overview", (
            "功能注释来源。", "最终基因中被 InterProScan 与 eggNOG-mapper 都注释、只被其一注释、或都没有注释的基因占比。"),
            values_text="{:.1f}%")
        fparts.append("<h3>总览</h3>" + table(oh, orows) + fig)
    h, r = read_tsv(O("06_Function/merged/function_stats.tsv"))
    if h:
        # Preferred_name / Description are eggNOG name and description fields written to the GFF3
        # (gene_name= / product=), not annotation categories: table only, not in the figure
        fields = {"Preferred_name": "基因名（gene_name，eggNOG Preferred_name）",
                  "Description": "功能描述（product，eggNOG Description）"}
        cats = [x for x in r if x[0] not in fields]
        rows_ = [[fields.get(x[0], x[0])] + x[1:] for x in r]
        fparts.append("<h3>各类注释的基因数</h3>" +
                      bar_figure([x[0] for x in cats], [float(x[2]) for x in cats], "% genes", figdir, "function_rates", (
                          "各类功能注释的覆盖率。", "有该类注释的基因占最终基因的百分比。any = 至少有一类；InterPro / Pfam 来自 InterProScan，"
                          "COG、KEGG、EC 来自 eggNOG-mapper，GO 为两者合并（灰色）。基因名和功能描述不是注释类别，只列在下表。"),
                          colors=[FUNC_SRC_COL.get(x[0], "#999999") for x in cats],
                          legend=[("InterProScan", "#2f6fb0"), ("eggNOG-mapper", "#e0a43a"), ("both (merged)", "#999999")],
                          values_text="{:.1f}%") +
                      table(h, rows_) +
                      "<p class='muted'>基因名和功能描述由 eggNOG-mapper 从直系同源基因转移而来，分别写入最终 GFF3 "
                      "（06_Function/merged/&lt;前缀&gt;.func.gff3）的 gene_name= 和 product= 属性。基因名只在同源基因有公认名称时才给出，"
                      "覆盖率低是正常的；功能描述在 eggNOG 找到同源组时一般都有，覆盖率通常与 COG 相同。</p>")
    h, r = read_tsv(F("interpro_databases"))
    if h:
        fparts.append("<h3>InterProScan 各成员数据库</h3>" +
                      bar_figure([x[0] for x in r if float(x[2]) > 0], [float(x[2]) for x in r if float(x[2]) > 0],
                                 "% genes", figdir, "interpro_databases", (
                                     "InterProScan 各成员数据库的命中率。", "被该数据库命中的基因占最终基因的百分比（没有命中的数据库不画）。"),
                                 values_text="{:.1f}%") +
                      table(h, r))
    for n_, title in (("top_pfam", "最常见的 Pfam 结构域"), ("top_interpro", "最常见的 InterPro 条目")):
        h, r = read_tsv(F(n_))
        if h and r:
            fparts.append(f"<h3>{title}</h3>" + table(h, r))
    gh, gr = read_tsv(F("go_namespace"))
    th, tr = read_tsv(F("go_top"))
    go_html = ""
    if gh and gr:
        go_html += table(gh, gr)
    elif th:
        go_html += "<p class='muted'>未提供 GO 本体文件（function.go_obo），只显示 GO 编号，不分 BP/MF/CC。</p>"
    if th and tr:
        go_html += table(th, tr)
    if go_html:
        fparts.append("<h3>GO</h3>" + go_html)
    kh, kr = read_tsv(F("kegg"))
    ph, pr = read_tsv(F("kegg_top_pathways"))
    if kh:
        fparts.append("<h3>KEGG</h3>" + table(kh, kr) +
                      ("<p class='muted'>最常见的通路（map 编号）：</p>" + table(ph, pr) if pr else ""))
    ch_, cr_ = read_tsv(F("cog_categories"))
    if ch_ and cr_:
        fparts.append("<h3>COG 功能分类（eggNOG-mapper）</h3>" +
                      bar_figure([f"{x[0]}: {x[1]}" for x in cr_], [int(x[2]) for x in cr_], "genes", figdir,
                                 "cog_categories", ("COG 功能分类。", "eggNOG-mapper 给出的 COG 类别的基因数；一个基因可属于多个类别，"
                                                    "S = 功能未知。"), wrap=48, values_text="{:.0f}") + table(ch_, cr_))
    lh, lr = read_tsv(F("eggnog_levels"))
    if lh and lr:
        fparts.append("<h3>eggNOG 注释所在的分类层级</h3>" + table(lh, lr))
    S.append(("function", "06_Function 功能注释", "".join(fparts)))
    # appendix
    vers = []
    for f in sorted(glob.glob(O("*/logs/versions/*.txt"))):
        txt = open(f, errors="replace").read().strip().splitlines()
        vers.append([os.path.relpath(f, a.outdir), (txt[0] if txt else "")[:120]])
    scripts = [[os.path.relpath(f, a.outdir)] for f in sorted(glob.glob(O("*/scripts/*.sh")))]
    cfg = html.escape(open(a.config).read()) if os.path.exists(a.config) else ""
    S.append(("appendix", "附录：软件版本、脚本索引、参数快照",
              "<h3>软件版本</h3>" + table(["文件", "版本信息（首行）"], vers) +
              "<h3>各步骤脚本</h3>" + table(["脚本"], scripts) +
              f"<h3>config</h3><pre>{cfg}</pre>"))

    banner = ""
    if a.no_rna:   # sections that need RNA-seq are not produced; say so prominently
        S = [x for x in S if x[0] not in ("rnaqc", "rnaalign", "eviann")]
        banner = ("<div class='note'><b>无 RNA 模式</b>（rna.allow_none: true）：未使用 RNA-seq。"
                  "BRAKER4 为 EP 模式（只用蛋白证据），不运行 EviAnn、StringTie，GeMoMa 不使用 RNA 剪接证据；"
                  "EVM 只整合 BRAKER4、GeMoMa 和 miniprot，权重调优只用蛋白证据打分；修正步骤没有 EviAnn 可比较，"
                  "所有位点按 R5 检查蛋白 / BRAKER4 证据；最终过滤以蛋白证据代替 RNA 证据。"
                  "结果没有 UTR、可变剪接和 lncRNA，基因结构的可靠性低于有 RNA 的注释。</div>")

    def nav_item(k, t):
        step, _, rest = t.partition(" ")
        if re.match(r"^\d\d_", step) and rest:                  # "05_EVM 修正（…）" -> step label + title
            return f"<a href='#{k}'><span class='step'>{html.escape(step)}</span>{html.escape(rest)}</a>"
        return f"<a href='#{k}'>{html.escape(t)}</a>"
    nav = ("<div class='navtitle'>目录</div>" + "".join(nav_item(k, t) for k, t, _ in S) +
           "<a class='top' href='#top'>↑ 回到顶部</a>")
    body = number_figures("".join(f"<section id='{k}'><h2>{html.escape(t)}</h2>{c}</section>" for k, t, c in S))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(f"<!doctype html><html lang='zh'><head><meta charset='utf-8'>"
                 f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
                 f"<title>{html.escape(a.prefix)} 注释报告</title><style>{CSS}</style></head>"
                 f"<body><div class='layout'><nav class='side'>{nav}</nav>"
                 f"<main><h1 id='top'>{html.escape(a.prefix)} 基因组注释统计报告</h1>"
                 f"<p class='muted'>由 Concordia 基因组注释流程（Snakemake）自动生成</p>{banner}{body}</main></div>"
                 f"<script>{NAV_JS}</script></body></html>")


if __name__ == "__main__":
    main()
