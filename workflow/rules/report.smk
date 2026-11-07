# 07_Report: illustrated statistics report covering every step.


# sources whose BUSCO / OMArk / AGAT appear in the EviAnn, BRAKER4 and homology sections
REPORT_QC = [k for k in QC_DIRS if k in ("eviann", "braker4", "gemoma") or k.startswith("gemoma_")]


rule report:
    input:
        rules.data_qc.input,
        O("01_Repeat/EDTA/genome.fa.mod.EDTA.TEanno.sum"),
        rules.repeat_stats.output,
        [O("02_RNA/qc/alignment_summary.tsv"), O("02_RNA/qc/merged_bam_summary.tsv"),
         O("02_RNA/EviAnn/stats/eviann_stats.tsv")] if HAS_RNA else [],
        O("03_Ab_initio/stats/braker4_stats.tsv"),
        O("03_Ab_initio/stats/braker4_support.tsv"),
        O("04_Homology/stats/gemoma_stats.tsv"),
        O("04_Homology/stats/gemoma_references.tsv"),
        O("04_Homology/stats/gemoma_attributes.tsv"),
        O("04_Homology/stats/miniprot_species.tsv"),
        rules.evm_stats.output,
        [f"{TUNE_DIR}/tuning_results.tsv"] if TUNE["enabled"] else [],
        rules.final_stats.output,
        O("05_EVM/final/filter_decisions.tsv"),
        rules.final_export.output,
        O("05_EVM/final/isoform_selection.tsv"),
        # per-source BUSCO / OMArk / AGAT shown in the EviAnn, BRAKER4 and GeMoMa sections
        [O(QC_DIRS[s], "agat/stats.txt") for s in REPORT_QC],
        [O(QC_DIRS[s], "omark/omark.sum") for s in REPORT_QC],
        [O(QC_DIRS[s], "busco", lin, "short_summary.txt") for s in REPORT_QC for lin in P.busco_lineages],
        O("05_EVM/qc/summary.tsv"),
        [O("05_EVM/qc/reference_comparison.tsv")] if P.reference_gff else [],
        O("05_EVM/refine/refine_decisions.tsv"),
        O("05_EVM/final/removed_genes.tsv"),
        O("06_Function/merged/function_stats.tsv"),
        rules.function_stats.output,
    output:
        O("07_Report/annotation_report.html"),
    params:
        script=script("07_Report", "01_report"),
        out=OUT,
        config=CONFIG_FILE,
        mode=P.braker4_mode(),
        tag=P.eviann_reads()[1],
        eviann_proteins=("未运行（无 RNA）" if not HAS_RNA else
                         "UniProt-SwissProt（eviann.swissprot）" if EVIANN_SWISSPROT
                         else "同源物种蛋白（未设置 eviann.swissprot）"),
        homologs=" ".join(HOMOLOG_NAMES),
        no_rna="" if HAS_RNA else "--no-rna",
    log:
        **logs("07_Report", "01_report"),
    conda:
        "../envs/core.yaml"
    shell:
        SH(PY("make_report.py") + """ --outdir {params.out} --prefix """ + PREFIX + """ \\
            --braker4-mode {params.mode} --eviann-tag {params.tag} --config {params.config} \
            --eviann-proteins '{params.eviann_proteins}' --homologs {params.homologs} {params.no_rna} --out {output}
        """, env="core")
