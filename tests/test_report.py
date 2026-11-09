"""07_Report must contain every genome assessment result that exists."""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "workflow", "scripts")


def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def build(out, with_optional=True):
    q = os.path.join(out, "00_Data", "qc")
    a = os.path.join(q, "assembly")
    w(f"{out}/00_Data/genome/rename_report.tsv", "renamed\tno\nreason\tok\n")
    w(f"{q}/seqkit/genome.stats.tsv", "file\tformat\tnum_seqs\nGENOME_X\tFASTA\t3\n")
    w(f"{a}/assembly_stats.tsv", "metric\tscaffold\tcontig\nN50\t123456\t65432\ngaps\t7\t-\n")
    w(f"{a}/sequences.tsv", "seq_id\tlength\tGC_percent\tN_bases\tgaps\tgap_length\tcontigs\n"
                            "SEQ_ALPHA\t200000\t35.0\t10\t1\t10\t2\nSEQ_BETA\t100000\t36.0\t0\t0\t0\t1\n")
    w(f"{q}/busco_genome/eudicots_odb12/short_summary.txt", "\tC:97.7%[S:95.1%,D:2.6%],F:0.8%,M:1.5%,n:2326\n")
    w(f"{a}/telomere/telomere_summary.tsv", "motif\tTTTAGGG\ntelomere_both_ends\t5\ntelomere_one_end\t2\n")
    w(f"{a}/telomere/telomere_per_sequence.tsv",
      "seq_id\tlength\tstart_repeats\tend_repeats\ttelomere_start\ttelomere_end\nSEQ_ALPHA\t200000\t310\t280\tyes\tyes\n")
    rows = ["section\tmetric\tvalue", "contiguity\tscaffold_N50\t123456",
            "completeness\tBUSCO_genome_eudicots_odb12\tC:97.7% [S:95.1%, D:2.6%], F:0.8%, M:1.5%, n:2326",
            "telomere\tmotif\tTTTAGGG"]
    if with_optional:
        w(f"{a}/merqury/genome.qv", "genome\t1000\t5000000\t43.21\t4.8e-05\n")
        w(f"{a}/merqury/genome.completeness.stats", "genome\tall\t900\t1000\t96.54\n")
        w(f"{a}/fcs_gx/fcs_gx_report.txt", "##[\"FCS genome report\"]\n#seq_id\tstart_pos\tend_pos\tseq_len\taction\n"
                                          "SEQ_CONTAM\t1\t5000\t5000\tEXCLUDE\n")
        rows += ["accuracy\tMerqury_QV\tQV=43.21", "contamination\tFCS-GX_flagged_regions\tEXCLUDE=1"]
    else:
        rows += ["accuracy\tMerqury_QV\tnot run", "contamination\tFCS-GX\tnot run"]
    w(f"{q}/assembly_quality.tsv", "\n".join(rows) + "\n")
    w(f"{out}/02_RNA/qc/alignment_summary.tsv",
      "sample\ttype\tplatform\treads\tmapped\tmapping_%\tunique_%\tmulti_%\tspliced_%\terror_rate\tmean_length\n"
      "LEAF_SR\tshort\tillumina\t20000\t19000\t95.00\t85.00\t5.00\t25.00\t2.1e-03\t150\n"
      "ISO_LR\tlong\tpacbio_isoseq\t1000\t980\t98.00\t-\t-\t90.00\t1.2e-02\t2500\n")
    w(f"{out}/02_RNA/qc/merged_bam_summary.tsv",
      "bam\ttotal_records\tprimary_mapped\tprimary_mapped_%\tused_by\nsr.bam\t20000\t17100\t90.00\tGeMoMa ERE\n")
    w(f"{out}/02_RNA/EviAnn/stats/eviann_stats.tsv",
      "metric\tcoding\tlncRNA\ngenes\t31415\t2718\ntranscripts_per_gene\t1.62\t1.05\n")
    w(f"{out}/02_RNA/EviAnn/stats/eviann_distributions.tsv",
      "set\tkind\tvalue\ncoding\tisoforms_per_gene\t2\ncoding\texons_per_transcript\t5\ncoding\tcds_aa\t300\n"
      "lncRNA\tisoforms_per_gene\t1\nlncRNA\texons_per_transcript\t2\n")
    w(f"{out}/02_RNA/EviAnn/qc/busco/eudicots_odb12/short_summary.txt", "\tC:93.3%[S:80.0%,D:13.3%],F:2.0%,M:4.7%,n:2326\n")
    w(f"{out}/02_RNA/EviAnn/qc/omark/omark.sum", "S:Single:90.1%,Duplicated:5.0%\nA:Consistent:88.8%\n")
    w(f"{out}/03_Ab_initio/stats/braker4_stats.tsv", "metric\tBRAKER4\ngenes\t27182\n")
    w(f"{out}/03_Ab_initio/stats/braker4_distributions.tsv", "set\tkind\tvalue\nBRAKER4\tcds_aa\t420\n")
    w(f"{out}/03_Ab_initio/stats/braker4_support.tsv",
      "column\tmodels\tsupported\tsupported_%\tmean\n# BRAKER4 mode dual; 27182 models\nRNA_SUPPORT_COL\t27182\t20000\t73.6\t0.6\n")
    w(f"{out}/03_Ab_initio/qc/busco/eudicots_odb12/short_summary.txt", "\tC:91.9%[S:89.0%,D:2.9%],F:3.0%,M:5.1%,n:2326\n")
    w(f"{out}/04_Homology/stats/gemoma_stats.tsv", "metric\tGeMoMa\ngenes\t24601\n")
    w(f"{out}/04_Homology/stats/gemoma_distributions.tsv", "set\tkind\tvalue\nGeMoMa\texons_per_transcript\t4\n")
    w(f"{out}/04_Homology/stats/gemoma_references.tsv",
      "species\tunfiltered_predictions\tspecies_final_transcripts\tmerged_representative\tmerged_supported\n"
      "REF_SPECIES_A\t30001\t20000\t15000\t18000\n")
    w(f"{out}/04_Homology/stats/gemoma_attributes.tsv", "metric\tvalue\ntie (intron RNA support) mean\t0.876\n")
    w(f"{out}/04_Homology/stats/miniprot_species.tsv",
      "species\tinput_proteins\taligned_proteins\taligned_%\talignments\tmean_identity\tidentity>=0.9\t"
      "with_frameshift_or_stop\nREF_SPECIES_A\t27000\t21000\t77.8\t25000\t0.654\t3000\t100\n")
    w(f"{out}/04_Homology/qc/GeMoMa_REF_SPECIES_A/busco/eudicots_odb12/short_summary.txt",
      "\tC:88.8%[S:86.0%,D:2.8%],F:3.0%,M:8.2%,n:2326\n")
    e = f"{out}/05_EVM"
    w(f"{e}/logs/evm_impl_check.txt", "standard\t/opt/bin/EVidenceModeler\tconfig evm.implementation=standard\n")
    w(f"{e}/inputs/weights.txt", "ABINITIO_PREDICTION\tBRAKER4\t2\nOTHER_PREDICTION\tEviAnn\t8\n")
    w(f"{e}/stats/evm_stats.tsv", "metric\tEVM\ngenes\t29001\n")
    w(f"{e}/stats/evm_distributions.tsv", "set\tkind\tvalue\nEVM\tcds_aa\t350\n")
    w(f"{e}/stats/evm_source_support.tsv", "source\ttranscripts\texact\texact_%\toverlap\toverlap_%\n"
                                          "EviAnn\t29001\t17777\t61.3\t25000\t86.2\n")
    w(f"{e}/stats/evm_source_combinations.tsv", "combination\ttranscripts\tpercent\nBRAKER4+GeMoMa+EviAnn\t12345\t42.6\n")
    w(f"{e}/final/filter_decisions.tsv",
      "# max_te_overlap=0.5 min_cds_aa=50; criterion c uses gene_support.tsv column(s): RNA\n"
      "gene\tsource\tcds_bp\tlongest_cds_aa\tte_overlap_bp\tte_overlap_frac\ttesorter\trna_stringtie\trna_eviann\t"
      "rna_braker4\trna_supported\taction\treason\n"
      "g1\tEVM\t600\t200\t480\t0.800\t-\tno\tno\tno\tno\tremove\tTE_overlap_0.80_no_RNA\n"
      "g2\tEVM\t600\t200\t480\t0.800\tLTR/Copia\tyes\tno\tno\tyes\tkeep\t-\n")
    w(f"{e}/final/removed_genes.tsv", "gene\treason\ng1\tTE_overlap_0.80_no_RNA\n")
    w(f"{e}/stats/final_stats.tsv", "metric\tfinal\tlncRNA\ngenes\t28123\t2718\n")
    w(f"{e}/stats/final_distributions.tsv", "set\tkind\tvalue\nfinal\tcds_aa\t360\n")
    w(f"{e}/stats/final_source_support.tsv", "source\ttranscripts\texact\texact_%\toverlap\toverlap_%\n"
                                            "GeMoMa\t28123\t15555\t55.3\t24000\t85.3\n")
    w(f"{e}/stats/final_source_combinations.tsv", "combination\ttranscripts\tpercent\nno exact match\t999\t3.6\n")
    w(f"{e}/stats/all_sets_stats.tsv", "metric\tBRAKER4\tGeMoMa\tEviAnn\tEVM\trefine\tfinal\n"
                                      "genes\t27182\t24601\t31415\t29001\t29500\t28123\n"
                                      "transcripts\t30000\t26000\t50000\t29001\t36000\t35000\n")
    w(f"{e}/qc/final/busco/eudicots_odb12/short_summary.txt", "\tC:98.4%[S:94.0%,D:4.4%],F:0.6%,M:1.0%,n:2326\n")
    w(f"{e}/final/AB.gff3", "c1\tx\tgene\t1\t10\t.\t+\t.\tID=AB1\n")
    w(f"{e}/final/AB.pep.fa", ">AB1.t1\nMA\n>AB2.t1\nMV\n")
    fs = f"{out}/06_Function/stats"
    w(f"{fs}/overview.tsv", "category\tgenes\tpercent\ntotal genes\t28123\t100.00\nboth\t20202\t71.83\n"
                           "InterProScan only\t1000\t3.56\neggNOG-mapper only\t900\t3.20\nneither (unannotated)\t6021\t21.41\n")
    w(f"{fs}/interpro_databases.tsv", "database\tgenes\tpercent\nPANTHER_DB_X\t19000\t67.56\n")
    w(f"{fs}/top_pfam.tsv", "pfam\tdescription\tgenes\nPF00069\tProtein kinase domain\t1111\n")
    w(f"{fs}/top_interpro.tsv", "interpro\tdescription\tgenes\nIPR011009\tKinase-like_dom_sf\t1200\n")
    w(f"{fs}/go_top.tsv", "go\tname\tnamespace\tgenes\nGO:0004672\tprotein kinase activity\tMF\t900\n")
    w(f"{fs}/go_namespace.tsv", "namespace\tgenes\tpercent\nBP\t14000\t49.78\n")
    w(f"{fs}/kegg.tsv", "category\tgenes\tpercent\nwith KEGG KO\t12121\t43.10\n")
    w(f"{fs}/kegg_top_pathways.tsv", "pathway\tname\tgenes\nmap01100\tMetabolic pathways\t2222\n")
    w(f"{fs}/cog_categories.tsv", "code\tdescription\tgenes\tpercent\nT\tSignal transduction mechanisms\t1717\t6.11\n")
    w(f"{fs}/eggnog_levels.tsv", "max_annot_lvl\tgenes\n33090|Viridiplantae\t15151\n")
    w(f"{out}/06_Function/merged/function_stats.tsv", "database\tgenes\tpercent\nany\t960\t85.41\n"
      "InterPro\t914\t81.32\nCOG_category\t918\t81.67\nPreferred_name\t63\t5.60\nDescription\t918\t81.67\n")
    rs = f"{out}/01_Repeat/stats"
    w(f"{rs}/overview.tsv", "item\tbp\tgenome_%\ngenome size (bp)\t10000\t\nTE annotated by EDTA (merged)\t4321\t43.21\n")
    w(f"{rs}/te_classes.tsv", "class\tsuperfamily\tcopies\tbp\tgenome_%\nLTR\t(all)\t30\t3000\t30.00\n"
                             "LTR\tGYPSY_SF\t20\t2000\t20.00\ntotal TE\t-\t30\t3000\t30.00\n")
    w(f"{rs}/te_methods.tsv", "method\tcopies\tbp\tgenome_%\nstructural\t7\t777\t7.77\n")
    w(f"{rs}/library.tsv", "classification\tconsensus_sequences\nLTR/Gypsy\t11\n")
    w(f"{rs}/landscape.tsv", "kimura_div\tclass\tbp\tgenome_%\n0\tLTR\t100\t1.00\n1\tDNA\t50\t0.50\n")
    w(f"{rs}/ltr_age.tsv", "# mutation rate u = 1.3e-8 (EDTA default, rice); T = K / (2u)\n"
                          "ltr\tsuperfamily\tidentity\tinsertion_mya\nc1:1..2\tGypsy\t0.99\t0.385\nc1:5..9\tGypsy\t0.98\t0.769\n")
    w(f"{out}/05_EVM/refine/refine_decisions.tsv",
      "# support=long min_support=2 min_cds_aa=50 max_te_overlap=0.5; criterion c not used\n"
      "locus\tseqid\tstart\tend\tstrand\tevm_genes\teviann_genes\trule\taction\tevm_support\teviann_support\tevidence\tnote\n"
      "L000001\tc1\t1\t900\t+\te1\tv1\tR1\tmerge_utr\t3\t3\t-\t-\n"
      "L000002\tc1\t1000\t2000\t+\te2\tvA,vB\tR4_split\tkeep_evm\t0\t0\tlinking_reads=0\tneeds_review: split not supported\n")
    w(f"{q}/homologs/ath.clean_report.tsv", "item\tcount\nkept genes\t27655\ndropped gene: type transposable_element_gene\t35090\n")
    w(f"{q}/homologs/osa.clean_report.tsv", "item\tcount\nkept genes\t35489\ndropped gene: type pseudogene\t2367\n")
    w(f"{e}/final/id_map.tsv", "feature\told_id\tnew_id\tseqid\tstart\tend\n"
                              "gene\tg1\tAB01G000010\tChr1\t1\t9\ngene\tg2\tAB01G000020\tChr1\t20\t30\n"
                              "mRNA\tg2.t\tAB01G000020.1\tChr1\t20\t30\ngene\tg3\tAB00G000010\tscaf15\t1\t9\n"
                              "gene\tn1\tAB01L000010\tChr1\t40\t50\n")
    w(f"{e}/tune/tuning_results.tsv", "# region Chr1 (Chr1, 50000 bp); P = multi-exon chains confirmed by any of lr_chain, prot_chain\n"
                                     "rank\tcombo\tweights\tgenes\tmulti_exon_tx\tlr_chain\tprot_chain\tlr_loci\tprot_loci\tP\tR\tF\n"
                                     "1\tc07\tBRAKER4=1, EviAnn=10\t500\t400\t0.8\t0.6\t0.9\t0.8\t0.9100\t0.8500\t0.8790\n"
                                     "2\tc00\tBRAKER4=2, EviAnn=8\t480\t390\t0.7\t0.6\t0.9\t0.8\t0.8800\t0.8500\t0.8647\n")
    w(f"{e}/tune/best_weights.txt", "ABINITIO_PREDICTION\tBRAKER4\t1\nOTHER_PREDICTION\tEviAnn\t10\n")
    w(f"{out}/cfg.yaml", "outdir: x\n")


def render_dir(tmp_path, out):
    html_ = str(tmp_path / "r.html")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_report.py"), "--outdir", out, "--prefix", "AB",
                        "--braker4-mode", "dual", "--eviann-tag", "bam_mix", "--config", f"{out}/cfg.yaml",
                        "--out-zh", html_], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return open(html_).read()


def render(tmp_path, **kw):
    out = str(tmp_path / "res")
    build(out, **kw)
    html = str(tmp_path / "r.html")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_report.py"), "--outdir", out, "--prefix", "AB",
                        "--braker4-mode", "dual", "--eviann-tag", "bam_mix", "--config", f"{out}/cfg.yaml",
                        "--eviann-proteins", "SWISSPROT_SET", "--homologs", "REF_SPECIES_A",
                        "--out-zh", html], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return open(html).read()


def test_all_assessment_results_in_report(tmp_path):
    h = render(tmp_path)
    for needle in ["基因组组装质量", "123456", "65432",           # summary + contiguity (scaffold/contig)
                   "SEQ_ALPHA", "GENOME_X",                     # per-sequence table, seqkit stats
                   "97.7", "eudicots_odb12",                    # BUSCO genome mode
                   "43.21", "96.54",                            # Merqury QV + completeness
                   "TTTAGGG", "telomere_both_ends",             # telomere summary + per-sequence
                   "SEQ_CONTAM", "EXCLUDE"]:                    # FCS-GX flagged region
        assert needle in h, needle
    assert "<h3>LAI</h3>" not in h                              # LAI off by default


def test_rna_alignment_in_report(tmp_path):
    h = render(tmp_path)
    for needle in ["RNA-seq 比对", "LEAF_SR", "ISO_LR", "95.00", "85.00", "pacbio_isoseq", "2500",
                   "sr.bam", "GeMoMa ERE", "rna_alignment"]:            # rna_alignment = figure alt text
        assert needle in h, needle


def test_not_run_items_are_listed(tmp_path):
    h = render(tmp_path, with_optional=False)
    assert "未运行" in h and "Merqury_QV" in h and "FCS-GX" in h


def test_eviann_stats_in_report(tmp_path):
    h = render(tmp_path)
    for needle in ["EviAnn 注释", "31415", "2718", "1.62", "bam_mix",   # stats table + reads tag
                   "93.3", "Consistent:88.8%", "eviann_distributions"]:  # BUSCO, OMArk, figure
        assert needle in h, needle


def test_braker4_and_gemoma_stats_in_report(tmp_path):
    h = render(tmp_path)
    for needle in ["BRAKER4 注释", "27182", "RNA_SUPPORT_COL", "73.6", "BRAKER4 mode dual", "91.9",
                   "braker4_distributions", "dual",
                   "同源注释（GeMoMa / miniprot）", "24601", "REF_SPECIES_A", "30001", "0.876", "gemoma_distributions",
                   "C:88.8%", "miniprot 蛋白比对", "0.654", "SWISSPROT_SET"]:
        assert needle in h, needle


def test_evm_and_final_stats_in_report(tmp_path):
    h = render(tmp_path)
    for needle in ["EVM 整合", "standard", "/opt/bin/EVidenceModeler", "BRAKER4</td><td>2", "29001", "17777",
                   "BRAKER4+GeMoMa+EviAnn", "evm_source_combinations",
                   "TE 过滤明细", "其中因有 RNA 证据而保留", "max_te_overlap=0.5",
                   "最终注释结果", "28123", "2718", "29500", "15555", "98.4", "final_distributions",
                   "05_EVM/final/AB.pep.fa", "2 序列", "1 基因", "基因结构对比"]:
        assert needle in h, needle


def test_function_stats_in_report(tmp_path):
    h = render(tmp_path)
    for needle in ["功能注释", "20202", "6021", "function_overview", "PANTHER_DB_X", "interpro_databases",
                   "PF00069", "Protein kinase domain", "IPR011009", "protein kinase activity", "14000",
                   "12121", "map01100", "Metabolic pathways", "Signal transduction mechanisms", "cog_categories", "Viridiplantae"]:
        assert needle in h, needle


def test_repeat_stats_in_report(tmp_path):
    h = render(tmp_path)
    for needle in ["01_Repeat 重复序列", "4321", "43.21", "GYPSY_SF", "te_classes", "777", "共 11 条一致序列",
                   "repeat_landscape", "1.3e-8 (EDTA default, rice)",
                   "<td>Gypsy</td><td>2</td>",
                   "只有 2 个完整 LTR"]:                 # too few intact LTRs: table only, no histogram
        assert needle in h, needle
    assert "ltr_insertion_time" not in h


def test_figures_have_numbered_captions(tmp_path):
    h = render(tmp_path)
    n_fig = h.count("<figure>")
    assert n_fig > 3 and h.count("<figcaption>") == n_fig
    assert "@@FIG@@" not in h
    assert [int(x) for x in re.findall(r"<b>图 (\d+)\. ", h)] == list(range(1, n_fig + 1))   # document order


def test_refine_section(tmp_path):
    h = render(tmp_path)
    for needle in ["R1 一致", "support=long min_support=2", "merge_utr", "需要人工检查的位点（1 个）", "L000002",
                   "refine_rules"]:
        assert needle in h, needle


def test_homolog_cleaning_in_report(tmp_path):
    h = render(tmp_path)
    for needle in ["同源物种注释预处理", "<th>ath</th>", "<th>osa</th>", "27655", "35090", "35489", "2367"]:
        assert needle in h, needle


def test_naming_table_in_report(tmp_path):
    h = render(tmp_path)
    assert "基因命名与各染色体基因数" in h
    assert "<td>01</td><td>Chr1</td><td>2</td><td>1</td><td>AB01G000010</td>" in h
    assert "<td>00</td><td>scaf15</td><td>1</td><td>0</td><td>AB00G000010</td>" in h


def test_evm_tuning_in_report(tmp_path):
    h = render(tmp_path)
    for needle in ["权重调优", "最优组合 <b>c07</b>（F = 0.8790）", "基线权重 c00 排第 2 位", "evm_tuning",
                   "全基因组 EVM 实际使用的权重", "<td>EviAnn</td><td>10</td>"]:
        assert needle in h, needle


def test_refine_rules_details_use_config_values(tmp_path):
    out = str(tmp_path / "res")
    build(out)
    with open(f"{out}/cfg.yaml", "w") as fh:
        fh.write("refine:\n  min_cds_aa: 80          # comment\n  remove_unsupported_evm: true\nnaming:\n  step: 10\n")
    html_ = str(tmp_path / "r.html")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_report.py"), "--outdir", out, "--prefix", "AB",
                        "--braker4-mode", "dual", "--eviann-tag", "bam_mix", "--config", f"{out}/cfg.yaml",
                        "--out-zh", html_], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    h = open(html_).read()
    assert "修正规则详解" in h and "min_cds_aa = 80" in h and "min_lr_support = 2" in h   # config value / default
    assert "remove_unsupported_evm = true" in h and "删除" in h
    assert "<div class='tbl'><table>" in h                                                # centred table wrapper


def test_function_rates_names_not_in_figure(tmp_path):
    h = render(tmp_path)
    i = h.index("<h3>各类注释的基因数</h3>")
    sec = h[i:h.index("<h3>", i + 10)]
    assert "function_rates" in sec and "基因名（gene_name，eggNOG Preferred_name）" in sec
    assert "功能描述（product，eggNOG Description）" in sec and "<td>Preferred_name</td>" not in sec


def test_side_navigation(tmp_path):
    h = render(tmp_path)
    assert "<nav class='side'>" in h and "<div class='layout'>" in h and "<script>" in h
    nav = h[h.index("<nav class='side'>"):h.index("</nav>")]
    ids = re.findall(r"<section id='([^']+)'>", h)
    assert ids and all(f"href='#{i}'" in nav for i in ids)                       # every section is listed
    assert "<span class='step'>01_Repeat</span>重复序列" in nav


def test_no_rna_banner_and_sections(tmp_path):
    out = str(tmp_path / "res")
    build(out)
    html_ = str(tmp_path / "r.html")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_report.py"), "--outdir", out, "--prefix", "AB",
                        "--braker4-mode", "EP", "--eviann-tag", "none", "--no-rna", "--out-zh", html_,
                        "--config", f"{out}/cfg.yaml"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    h = open(html_).read()
    assert "无 RNA 模式" in h and "没有 UTR、可变剪接和 lncRNA" in h
    for sec in ("rnaqc", "rnaalign", "eviann"):
        assert f"<section id='{sec}'>" not in h and f"href='#{sec}'" not in h
    assert "<section id='braker4'>" in h


def test_long_read_weights_warning(tmp_path):
    out = str(tmp_path / "res")
    build(out)
    os.makedirs(f"{out}/02_RNA/qc", exist_ok=True)
    with open(f"{out}/02_RNA/qc/long_read_weights.tsv", "w") as fh:
        fh.write("sample\tplatform\taligned_reads\tweighted_reads\tweighted_%\tfull_length_reads\tweight_source\n"
                 "isoA\tpacbio_isoseq\t100\t0\t0.0\t100\tfull_length_coverage in read headers\n"
                 "isoB\tpacbio_isoseq\t100\t100\t100.0\t900\tcount table c.csv\n")
    h = render_dir(tmp_path, out)
    assert "长读的全长 reads 权重" in h and "Iso-Seq 样本 isoA 没有任何一条记录带全长 reads 数" in h
    assert "isoB 没有" not in h


def render_both(tmp_path):
    out = str(tmp_path / "res")
    build(out)
    zh, en = str(tmp_path / "annotation_report_zh.html"), str(tmp_path / "annotation_report_en.html")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_report.py"), "--outdir", out, "--prefix", "AB",
                        "--braker4-mode", "dual", "--eviann-tag", "6 x bam + 1 x bam_isoseq", "--config",
                        f"{out}/cfg.yaml", "--eviann-proteins", "swissprot", "--homologs", "REF_SPECIES_A",
                        "--out-zh", zh, "--out-en", en], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return open(zh).read(), open(en).read()


def test_english_report(tmp_path):
    zh, en = render_both(tmp_path)
    assert "<html lang='en'>" in en and "<html lang='zh'>" in zh
    # same data, same sections and figures in both versions
    assert re.findall(r"<section id='([^']+)'>", zh) == re.findall(r"<section id='([^']+)'>", en)
    assert zh.count("<figure>") == en.count("<figure>") > 3
    assert [int(x) for x in re.findall(r"<b>Figure (\d+)\. ", en)] == list(range(1, en.count("<figure>") + 1))
    for needle in ["Genome assembly quality", "Refinement rules in detail", "Weight tuning", "Best set <b>c07</b>",
                   "UniProt-SwissProt (eviann.swissprot)", "6 x bam + 1 x bam_isoseq", "Contents"]:
        assert needle in en, needle
    assert "UniProt-SwissProt（eviann.swissprot）" in zh
    # links to each other
    assert "href='annotation_report_zh.html'" in en and "href='annotation_report_en.html'" in zh
    # no Chinese left in the English version (apart from the link to the Chinese one; data tables excluded)
    text = re.sub(r"<pre>.*?</pre>|src='data:[^']*'", "", en, flags=re.S).replace("中文版", "")
    assert not re.search(r"[一-鿿　-〿＀-￯]", re.sub(r"<[^>]+>", " ", text))


def test_runtime_and_resources(tmp_path):
    out = str(tmp_path / "res")
    build(out)
    bm = "s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n"
    w(f"{out}/01_Repeat/logs/benchmarks/01_edta.tsv", bm + "7200.5\t2:00:00\t30000.0\t1\t1\t1\t0\t0\t300\t28800.0\n")
    w(f"{out}/02_RNA/logs/benchmarks/02_hisat2.leaf.tsv", bm + "600\t0:10:00\t3000\t1\t1\t1\t0\t0\t50\t1200\n")
    w(f"{out}/02_RNA/logs/benchmarks/02_hisat2.root.tsv", bm + "900\t0:15:00\t3500\t1\t1\t1\t0\t0\t50\t1800\n")
    req = f"{out}/07_Report/resource_requests.tsv"
    w(req, "rule\tthreads\tmem_mb\truntime_min\tbenchmark\n"
           "edta\t48\t32768\t120\t01_Repeat/logs/benchmarks/01_edta.tsv\n"
           "hisat2_align\t16\t32000\t\t02_RNA/logs/benchmarks/02_hisat2.{sample,[^/]+}.tsv\n")
    en, zh = str(tmp_path / "en.html"), str(tmp_path / "zh.html")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_report.py"), "--outdir", out, "--prefix", "AB",
                        "--braker4-mode", "dual", "--eviann-tag", "x", "--config", f"{out}/cfg.yaml",
                        "--requests", req, "--out-zh", zh, "--out-en", en], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    h = open(en).read()
    sec = h[h.index("<section id='resources'>"):]
    assert "3 jobs" in sec and "runtime_top" in sec
    # edta: 4 cores of 48 used, 29.3 GB of 32 GB (> 90%), 2 h of 2 h requested
    assert "<td>edta</td><td>1</td><td>2:00:00</td>" in sec and "<td>48</td><td>4.0</td><td>8</td>" in sec
    assert "memory near or above the request: use ≥ 39G" in sec and "fewer threads" in sec and "raise runtime" in sec
    # hisat2: two jobs matched through the wildcard pattern; memory far below the request
    assert "<td>hisat2_align</td><td>2</td><td>0:15:00</td><td>0:25:00</td>" in sec
    assert "about 6G is enough" in sec
    assert "运行时间与资源" in open(zh).read()
