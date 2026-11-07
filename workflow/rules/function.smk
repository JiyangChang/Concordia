# 06_Function: InterProScan (chunked) + eggNOG-mapper on the final longest-isoform proteins.

FN = config.get("function", {}) or {}
N_CHUNKS = int(FN.get("chunks", 20))
CHUNKS = [f"{i:03d}" for i in range(1, N_CHUNKS + 1)]


rule func_split:
    input:
        FINAL_LONGEST_PEP,
    output:
        # temporary: deleted once every InterProScan chunk has run
        temp(expand(O("06_Function/interproscan/chunks/chunk_{chunk}.faa"), chunk=CHUNKS)),
    params:
        script=script("06_Function", "01_split"),
    log:
        **logs("06_Function", "01_split"),
    shell:
        SH(PY("split_fasta.py") + """ {input} {output}
        """)


IPS_DIR = P.path(FN.get("interproscan_dir")) or ""


rule interproscan_setup:
    """One-time InterProScan setup after unpacking the release (indexes the HMM models with
    setup.py, as required by the InterProScan 5 installation instructions). Skipped when done."""
    input:
        # run only when InterProScan is actually needed (proteins split), not at workflow start
        rules.func_split.output,
    output:
        O("06_Function/interproscan/setup.done"),
    params:
        script=script("06_Function", "01b_interproscan_setup"),
        ips=IPS_DIR,
    threads: 4
    log:
        **logs("06_Function", "01b_interproscan_setup"),
    conda:
        "../envs/function.yaml"
    shell:
        SH("""
        if [ ! -x {params.ips}/interproscan.sh ]; then echo "ERROR: {params.ips}/interproscan.sh not found" >&2; exit 1; fi
        cd {params.ips}
        if [ ! -f .pipeline_setup_done ]; then
            python3 setup.py -f interproscan.properties
            touch .pipeline_setup_done
        fi
        java -version 2>&1 | head -n 1
        touch {output}
        """, env="function")


rule interproscan:
    input:
        setup=rules.interproscan_setup.output,
        faa=O("06_Function/interproscan/chunks/chunk_{chunk}.faa"),
    output:
        # both declared (so a re-run deletes old files first; InterProScan would otherwise write chunk_*_1.gff3)
        # and temporary: deleted once interproscan_merge has combined them
        tsv=temp(O("06_Function/interproscan/chunk_{chunk}.tsv")),
        gff=temp(O("06_Function/interproscan/chunk_{chunk}.gff3")),
    params:
        script=script("06_Function", "02_interproscan.{chunk}"),
        exe=os.path.join(IPS_DIR, "interproscan.sh"),
        appl=f"--appl {FN['interproscan_appl']}" if FN.get("interproscan_appl") else "",
        base=O("06_Function/interproscan/chunk_{chunk}"),
        tmp=O("06_Function/interproscan/tmp_{chunk}"),
    threads: res("interproscan")["threads"]
    resources:
        mem_mb=res("interproscan")["mem_mb"],
        runtime=res("interproscan")["runtime"],
    log:
        **logs("06_Function", "02_interproscan.{chunk}"),
    conda:
        "../envs/function.yaml"
    shell:
        SH("""
        if [ ! -s {input.faa} ]; then : > {output.tsv}; : > {output.gff}; exit 0; fi
        rm -f {params.base}.tsv {params.base}.gff3 {params.base}_*.gff3 {params.base}_*.tsv
        mkdir -p {params.tmp}
        {params.exe} -i {input.faa} -f tsv,gff3 -goterms -pa -cpu {threads} -dp {params.appl} \\
                     -b {params.base} -T {params.tmp}
        rm -rf {params.tmp}
        """, env="function", versions=f"{IPS_DIR}/interproscan.sh -version 2>/dev/null | head -n 3")   # not on PATH: use the configured release


IPS_TSV = O("06_Function/interproscan", f"{PREFIX}.interproscan.tsv")
IPS_GFF = O("06_Function/interproscan", f"{PREFIX}.interproscan.gff3")


rule interproscan_merge:
    """Complete InterProScan result (all member databases, positions, E-values) as one TSV and one GFF3;
    the per-chunk files are temporary and deleted afterwards."""
    input:
        tsv=expand(O("06_Function/interproscan/chunk_{chunk}.tsv"), chunk=CHUNKS),
        gff=expand(O("06_Function/interproscan/chunk_{chunk}.gff3"), chunk=CHUNKS),
    output:
        tsv=IPS_TSV,
        gff=IPS_GFF,
    params:
        script=script("06_Function", "02b_interproscan_merge"),
    log:
        **logs("06_Function", "02b_interproscan_merge"),
    shell:
        SH(PY("merge_interproscan.py") + """ --tsv {input.tsv} --gff {input.gff} --out-tsv {output.tsv} --out-gff {output.gff}
        """)


rule eggnog:
    input:
        FINAL_LONGEST_PEP,
    output:
        O("06_Function/eggnog", f"{PREFIX}.emapper.annotations"),
    params:
        script=script("06_Function", "03_eggnog"),
        outdir=O("06_Function/eggnog"),
        db=P.path(FN.get("eggnog_data_dir")) or "",
        scope=FN.get("eggnog_tax_scope", "auto"),
    threads: res("eggnog")["threads"]
    resources:
        mem_mb=res("eggnog")["mem_mb"],
        runtime=res("eggnog")["runtime"],
    log:
        **logs("06_Function", "03_eggnog"),
    conda:
        "../envs/function.yaml"
    shell:
        SH("""
        emapper.py -i {input} --itype proteins -m diamond --data_dir {params.db} --cpu {threads} \\
                   --tax_scope {params.scope} -o """ + PREFIX + """ --output_dir {params.outdir} --override
        """, env="function", versions="emapper.py --version")


rule merge_function:
    input:
        gff=FINAL_GFF,
        ipr=IPS_TSV,
        egg=rules.eggnog.output,
    output:
        gff=O("06_Function/merged", f"{PREFIX}.func.gff3"),
        tsv=O("06_Function/merged/functional_annotation.tsv"),
        go=O("06_Function/merged/go_annotation.tsv"),
        kegg=O("06_Function/merged/kegg_annotation.tsv"),
        stats=O("06_Function/merged/function_stats.tsv"),
    params:
        script=script("06_Function", "04_merge"),
    log:
        **logs("06_Function", "04_merge"),
    shell:
        SH(PY("merge_function.py") + """ --gff {input.gff} --eggnog {input.egg} --interproscan {input.ipr} \\
            --out-gff {output.gff} --out-tsv {output.tsv} --out-go {output.go} --out-kegg {output.kegg} \\
            --out-stats {output.stats}
        cat {output.stats}
        """)


FUNC_STATS = ["overview", "interpro_databases", "top_pfam", "top_interpro", "go_top", "go_namespace",
              "kegg", "kegg_top_pathways", "cog_categories", "eggnog_levels"]


KEGG_NAMES_CFG = FN.get("kegg_pathway_names", "download")


rule kegg_pathway_names:
    """KEGG pathway names (map ID -> name): downloaded from the KEGG REST API, or a user file.
    A failed download leaves an empty table (the report then shows map IDs only)."""
    output:
        O("06_Function/stats/kegg_pathway_names.tsv"),
    params:
        script=script("06_Function", "05a_kegg_pathway_names"),
        src="" if KEGG_NAMES_CFG in (None, "", "download") else P.path(KEGG_NAMES_CFG),
        download="yes" if KEGG_NAMES_CFG == "download" else "no",
    log:
        **logs("06_Function", "05a_kegg_pathway_names"),
    shell:
        SH("""
        if [ -n "{params.src}" ]; then
            cp {params.src} {output}
        elif [ "{params.download}" = "yes" ] && curl -sSf --max-time 120 https://rest.kegg.jp/list/pathway -o {output}.tmp; then
            mv {output}.tmp {output}
        else
            rm -f {output}.tmp
            : > {output}
            echo "WARNING: no KEGG pathway names (download disabled or failed); report shows map IDs only" >&2
        fi
        echo "pathway names: $(wc -l < {output})"
        """)


rule function_stats:
    """Gene-level statistics of the functional annotation, for 07_Report."""
    input:
        gff=FINAL_GFF,
        ipr=IPS_TSV,
        egg=rules.eggnog.output,
        kegg_names=rules.kegg_pathway_names.output,
    output:
        [O(f"06_Function/stats/{n}.tsv") for n in FUNC_STATS],
    params:
        script=script("06_Function", "05_function_stats"),
        outdir=O("06_Function/stats"),
        obo=f"--go-obo {P.path(FN['go_obo'])}" if P.path(FN.get("go_obo")) else "",
        top=FN.get("report_top", 20),
    log:
        **logs("06_Function", "05_function_stats"),
    shell:
        SH(PY("function_stats.py") + """ --gff {input.gff} --interproscan {input.ipr} --eggnog {input.egg} \\
            {params.obo} --kegg-names {input.kegg_names} --top {params.top} --outdir {params.outdir}
        cat {params.outdir}/overview.tsv
        """)
