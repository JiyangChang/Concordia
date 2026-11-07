# 03_Ab_initio: BRAKER4 (nested Snakemake). Only PacBio long reads are passed to BRAKER4.

B4 = config.get("braker4", {}) or {}
B4_REPO = project_path(B4.get("repo", "external/BRAKER4"))
B4_DIR = O("03_Ab_initio/BRAKER4")
B4_BAM, B4_ISO = P.braker4_bams()
B4_MODE = P.braker4_mode()
B4_RESULT = O("03_Ab_initio/BRAKER4/output", PREFIX, "results")
# Snakemake of the running workflow (absolute interpreter path, also valid on SLURM nodes with a shared file system)
SNAKEMAKE = f"{sys.executable} -m snakemake"


rule braker4_inputs:
    input:
        genome=GENOME,
        masked=SOFTMASKED,
        proteins=PROTEINS,
        bams=[b for b in (B4_BAM, B4_ISO) if b],
    output:
        csv=O("03_Ab_initio/BRAKER4/samples.csv"),
        ini=O("03_Ab_initio/BRAKER4/config.ini"),
    params:
        script=script("03_Ab_initio", "01_braker4_inputs"),
        bam=B4_BAM,
        iso=B4_ISO,
        lineage=P.busco_lineages[0],
        # BUSCO datasets for BRAKER4's compleasm: qc.busco_download_path (shared with our BUSCO), else a directory in
        # the results (never the BRAKER4 code directory)
        busco_dl=P.busco_download_path or O("03_Ab_initio/BRAKER4/busco_downloads"),
        # the repository ships config.ini.example (no config.ini)
        template=next((os.path.join(B4_REPO, f) for f in ("config.ini", "config.ini.example")
                       if os.path.exists(os.path.join(B4_REPO, f))), os.path.join(B4_REPO, "config.ini.example")),
        min_contig=B4.get("min_contig", 10000),
        skip_busco=1 if B4.get("skip_busco", True) else 0,
        run_omark=1 if B4.get("run_omark", False) else 0,
        minisplice=1 if B4.get("use_minisplice", False) else 0,
        threads=res("braker4")["threads"],
        mem=res("braker4")["mem_mb"],
        runtime=res("braker4")["runtime"],
    log:
        **logs("03_Ab_initio", "01_braker4_inputs"),
    shell:
        SH("# BRAKER4 mode decided from the data: " + B4_MODE + "\n" + PY("make_braker4_inputs.py") + """ \\
            --sample """ + PREFIX + """ --genome {input.genome} --masked {input.masked} --proteins {input.proteins} \\
            --bam '{params.bam}' --isoseq-bam '{params.iso}' --busco-lineage {params.lineage} \\
            --busco-download-path {params.busco_dl} \\
            --min-contig {params.min_contig} --skip-busco {params.skip_busco} --run-omark {params.run_omark} \\
            --use-minisplice {params.minisplice} --threads {params.threads} --mem-mb {params.mem} \\
            --runtime {params.runtime} --template {params.template} \\
            --out-csv {output.csv} --out-ini {output.ini}
        """)


def braker4_executor_args():
    if B4.get("executor", "local") != "slurm":
        return ""
    s = B4.get("slurm", {}) or {}
    dr = " ".join(f"{k}={v}" for k, v in (s.get("default_resources") or {}).items())
    return f"--executor slurm --jobs {B4.get('max_jobs', 50)}" + (f" --default-resources {dr}" if dr else "")


rule braker4:
    input:
        csv=rules.braker4_inputs.output.csv,
        ini=rules.braker4_inputs.output.ini,
    output:
        gff=os.path.join(B4_RESULT, "braker.gff3.gz"),
    params:
        script=script("03_Ab_initio", "02_braker4"),
        workdir=B4_DIR,
        snakefile=os.path.join(B4_REPO, "Snakefile"),
        snakemake=SNAKEMAKE,
        cache=os.path.abspath(B4.get("container_cache", ".singularity_cache")),
        # 00_Data may hold symbolic links to the inputs: their directories must be visible in the containers;
        # so must the shared BUSCO dataset directory used by BRAKER4's compleasm
        bind=",".join(dict.fromkeys([x for x in str(B4.get("bind", "")).split(",") if x] + [OUT] +
                                    [os.path.dirname(os.path.realpath(f)) for f in (P.genome, P.orthodb) if f] +
                                    ([os.path.realpath(P.busco_download_path)] if P.busco_download_path else []))),
        executor=braker4_executor_args(),
    threads: res("braker4")["threads"]
    resources:
        mem_mb=res("braker4")["mem_mb"],
        runtime=res("braker4")["runtime"],
    log:
        **logs("03_Ab_initio", "02_braker4"),
    shell:
        # The nested BRAKER4 run uses the Snakemake that runs this workflow (same version 8.18.2 and SLURM plugin);
        # its tools come from BRAKER4's own containers, so no conda environment is needed here.
        SH("""
        export BRAKER4_CONFIG={input.ini}
        {params.snakemake} --snakefile {params.snakefile} --cores {threads} {params.executor} \\
                  --use-singularity --singularity-prefix {params.cache} --singularity-args "-B {params.bind}" \\
                  --latency-wait 120 --restart-times 2 --rerun-incomplete
        """, workdir="{params.workdir}", versions=f"{SNAKEMAKE} --version")


rule braker4_gff:
    input:
        rules.braker4.output.gff,
    output:
        gff=QC_GFF["braker4"],
        support=O("03_Ab_initio/braker4.gene_support.tsv"),  # empty if BRAKER4 did not write one
    params:
        script=script("03_Ab_initio", "03_braker4_gff"),
        support=os.path.join(B4_RESULT, "gene_support.tsv"),
    log:
        **logs("03_Ab_initio", "03_braker4_gff"),
    shell:
        SH("""
        gzip -dc {input} > {output.gff}
        if [ -f {params.support} ]; then cp {params.support} {output.support}; else : > {output.support}; fi
        """)


rule braker4_stats:
    """Gene-model statistics of the BRAKER4 annotation + gene_support.tsv summary, for 07_Report."""
    input:
        gff=rules.braker4_gff.output.gff,
        support=rules.braker4_gff.output.support,
    output:
        stats=O("03_Ab_initio/stats/braker4_stats.tsv"),
        dist=O("03_Ab_initio/stats/braker4_distributions.tsv"),
        support=O("03_Ab_initio/stats/braker4_support.tsv"),
    params:
        script=script("03_Ab_initio", "04_braker4_stats"),
        mode=B4_MODE,
    log:
        **logs("03_Ab_initio", "04_braker4_stats"),
    shell:
        SH(PY("annotation_stats.py") + """ --set BRAKER4={input.gff} --out {output.stats} --out-dist {output.dist}
        """ + PY("braker4_support.py") + """ --support {input.support} --mode {params.mode} --out {output.support}
        cat {output.stats} {output.support}
        """)
