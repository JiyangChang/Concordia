# 03_Ab_initio: BRAKER4 (nested Snakemake). Short reads (sr.bam) and PacBio long reads (pacbio.bam) are passed to
# BRAKER4 (mode dual / ETP / IsoSeq; EP with ONT only or without RNA); ONT long reads are not.

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
        # as busco_ready: left out only for BRAKER4 runs that exist without the flag (older analyses)
        busco=lambda wc: ancient(O("00_Data/busco_datasets/compleasm.ready"))
        if os.path.exists(O("00_Data/busco_datasets/compleasm.ready"))
        or not os.path.exists(O("03_Ab_initio/BRAKER4/samples.csv")) else [],
    output:
        csv=O("03_Ab_initio/BRAKER4/samples.csv"),
        ini=O("03_Ab_initio/BRAKER4/config.ini"),
    params:
        script=script("03_Ab_initio", "01_braker4_inputs"),
        bam=B4_BAM,
        iso=B4_ISO,
        lineage=P.busco_lineages[0],
        # BUSCO datasets for BRAKER4's compleasm, prepared by busco_compleasm_data: qc.busco_download_path when
        # writable (shared with our BUSCO), else 00_Data/busco_downloads in the results
        busco_dl=COMPLEASM_DIR,
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


# braker4.bind: comma-separated mount points, each "path" or "path:container path[:opts]"; host paths resolved like
# every other path setting ('~', relative to the working directory)
B4_BIND = [":".join([P.path(b.split(":")[0])] + b.split(":")[1:])
           for b in (x.strip() for x in str(B4.get("bind") or "").split(",")) if b]
# Apptainer refuses to start a container when a mounted path does not exist: warn at start-up (dry-runs keep
# working, e.g. while a disk is not mounted); the braker4 job checks again and stops with a clear message
_b4_missing = [b.split(":")[0] for b in B4_BIND if not os.path.exists(b.split(":")[0])]
if _b4_missing:
    sys.stderr.write("WARNING: braker4.bind: path(s) do not exist: " + ", ".join(_b4_missing) +
                     " - BRAKER4 will fail unless they exist when it runs (set existing directories or \"\")\n")


def braker4_executor_args():
    if B4.get("executor", "local") != "slurm":
        return ""
    s = B4.get("slurm", {}) or {}
    # --default-resources values are evaluated as Python expressions: text (partition, account ...) is quoted so that
    # names such as "all" stay text; numbers are passed as they are; already quoted values are kept
    def item(k, v):
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            return f"{k}={v}"
        v = str(v).strip()
        if re.fullmatch(r"[0-9.eE+-]+", v):
            return f"{k}={v}"
        if len(v) >= 2 and v[0] == v[-1] and v[0] in "'\"":
            v = v[1:-1]
        return shlex.quote(f"{k}={v!r}")
    dr = " ".join(item(k, v) for k, v in (s.get("default_resources") or {}).items())
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
        cache=P.path(B4.get("container_cache") or ".singularity_cache"),
        # added automatically: the results; the directories of the inputs 00_Data links to; the BUSCO dataset
        # directory of BRAKER4's compleasm; the BRAKER4 code (its scripts/ are run inside the containers)
        bind=",".join(dict.fromkeys(B4_BIND + [OUT] +
                                    [os.path.dirname(os.path.realpath(f)) for f in (P.genome, P.orthodb) if f] +
                                    [os.path.realpath(COMPLEASM_DIR), os.path.realpath(B4_REPO)])),
        bind_user=" ".join(b.split(":")[0] for b in B4_BIND),
        # BUSCO data of compleasm checked again before BRAKER4 starts (fetched if it went missing since the flag)
        busco_check=PY("fetch_busco.py") + f" --dir {COMPLEASM_DIR} --compleasm"
                    + (f" --lineage {P.busco_lineages[0]}" if P.busco_lineages else ""),
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
        for d in {params.bind_user}; do
            [ -e "$d" ] || {{ echo "ERROR: braker4.bind: $d does not exist (Apptainer cannot mount it): list existing directories, or leave it empty" >&2; exit 1; }}
        done
        {params.busco_check}
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
