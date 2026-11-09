# 00_Data: collect inputs, rename sequence IDs only when needed, assess raw data.

RN = config.get("rename_seqids", {}) or {}


rule genome_input:
    input:
        P.genome,
    output:
        GENOME_RAW,
    params:
        script=script("00_Data", "01_genome_input"),
    log:
        **logs("00_Data", "01_genome_input"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        # link the input when it is an uncompressed FASTA with regular line widths (what faidx needs);
        # otherwise (gzip, irregular wrapping) write a reformatted copy
        if gzip -t {input} 2> /dev/null; then
            echo "gzip input: decompressing"
            seqkit seq -w 60 {input} > {output}
        elif awk '/^>/ {{if (n > 1 && bad) exit 1; w = 0; bad = 0; n = 0; last = 0; next}}
                  {{n++; if (last || length($0) == 0) bad = 1; if (!w) w = length($0); else if (length($0) > w) bad = 1; else if (length($0) < w) last = 1}}
                  END {{if (n > 1 && bad) exit 1}}' {input}; then
            echo "plain FASTA with regular line widths: symbolic link"
            ln -sfn {input} {output}
        else
            echo "irregular line widths: reformatting"
            seqkit seq -w 60 {input} > {output}
        fi
        """, env="core", versions="seqkit version")


rule rename_seqids:
    input:
        GENOME_RAW,
    output:
        fasta=GENOME,
        map=SEQID_MAP,
        report=O("00_Data/genome/rename_report.tsv"),
    params:
        script=script("00_Data", "02_rename_seqids"),
        max_len=RN.get("max_len", 13),
        prefix=RN.get("prefix", "seq"),
        invalid="true" if RN.get("on_invalid_chars", True) else "false",
        min_len=config.get("min_contig_len", 0),
    log:
        **logs("00_Data", "02_rename_seqids"),
    shell:
        SH(PY("rename_seqids.py") + """ \\
            --fasta {input} --out-fasta {output.fasta} --out-map {output.map} --out-report {output.report} \\
            --max-len {params.max_len} --prefix {params.prefix} --min-len {params.min_len} \\
            --on-invalid-chars {params.invalid} --link-unchanged
        cat {output.report}
        """)


rule genome_stats:
    input:
        GENOME,
    output:
        O("00_Data/qc/seqkit/genome.stats.tsv"),
    params:
        script=script("00_Data", "03_genome_stats"),
    log:
        **logs("00_Data", "03_genome_stats"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        seqkit stats -a -T {input} > {output}
        """, env="core")


# BUSCO datasets: present ones are used, missing ones are downloaded first by busco_dataset (see lib/busco_data.py
# for where); BUSCO itself then always runs offline.
BUSCO_LOCAL = O("00_Data/busco_downloads")
COMPLEASM_DIR = busco_data.compleasm_dir(P.busco_download_path, BUSCO_LOCAL)
_busco_errs = busco_data.check_names(P.busco_lineages, [P.busco_download_path, BUSCO_LOCAL])
if _busco_errs:
    raise ConfigError("Config errors:\n  - " + "\n  - ".join(_busco_errs))


def busco_dir(lineage):
    return busco_data.lineage_dir(P.busco_download_path, BUSCO_LOCAL, lineage)


def busco_ready(output_pattern):
    """Input function: the dataset flag (ancient: a re-created flag triggers nothing). Left out only when the BUSCO
    result exists without a flag (analyses started before this step existed), so that such results are not re-run;
    otherwise always given, so that the set of inputs stays the same from run to run."""
    def f(wildcards):
        flag = O("00_Data/busco_datasets", f"{wildcards.lineage}.busco.ready")
        if os.path.exists(flag) or not os.path.exists(output_pattern.format(**wildcards)):
            return ancient(flag)
        return []
    return f


def busco_check(wildcards):
    """Run at the start of every BUSCO job: the .ready flag only says the dataset was there once; this checks it is
    still present (and complete) and fetches it again otherwise. No network access when it is there."""
    return PY("fetch_busco.py") + f" --dir {busco_dir(wildcards.lineage)} --lineage {wildcards.lineage}"


def busco_offline(wildcards):
    return f"--offline --download_path {busco_dir(wildcards.lineage)}"


localrules: busco_dataset, busco_compleasm_data, busco_datasets


rule busco_dataset:
    """Check that a BUSCO dataset is present; download it when missing (runs in the main Snakemake process, i.e. on
    the machine / login node that started the workflow, which usually has internet access)."""
    output:
        O("00_Data/busco_datasets/{lineage}.busco.ready"),
    params:
        script=script("00_Data", "00_busco_dataset.{lineage}"),
        dir=lambda wc: busco_dir(wc.lineage),
    log:
        **logs("00_Data", "00_busco_dataset.{lineage}"),
    shell:
        SH(PY("fetch_busco.py") + """ --dir {params.dir} --lineage {wildcards.lineage} --done {output}
        """)


rule busco_compleasm_data:
    """Files BRAKER4's compleasm needs (first BUSCO lineage, eukaryota_odb12, placement files, markers)."""
    output:
        O("00_Data/busco_datasets/compleasm.ready"),
    params:
        script=script("00_Data", "00_busco_compleasm"),
        dir=COMPLEASM_DIR,
        lineage=P.busco_lineages[0] if P.busco_lineages else "",
    log:
        **logs("00_Data", "00_busco_compleasm"),
    shell:
        SH(PY("fetch_busco.py") + """ --dir {params.dir} --lineage {params.lineage} --compleasm --done {output}
        """)


rule busco_datasets:
    """Target to fetch all BUSCO datasets in advance (e.g. on a login node): snakemake ... busco_datasets"""
    input:
        expand(O("00_Data/busco_datasets/{lineage}.busco.ready"), lineage=P.busco_lineages),
        O("00_Data/busco_datasets/compleasm.ready"),


rule busco_genome:
    input:
        genome=GENOME,
        ready=busco_ready(O("00_Data/qc/busco_genome/{lineage}/short_summary.txt")),
    output:
        O("00_Data/qc/busco_genome/{lineage}/short_summary.txt"),
    params:
        script=script("00_Data", "04_busco_genome.{lineage}"),
        outpath=O("00_Data/qc/busco_genome"),
        offline=busco_offline,
        check=busco_check,
    threads: res("busco")["threads"]
    resources:
        mem_mb=res("busco")["mem_mb"],
        runtime=res("busco")["runtime"],
    log:
        **logs("00_Data", "04_busco_genome.{lineage}"),
    conda:
        "../envs/qc.yaml"
    shell:
        SH("""
        {params.check}
        busco -i {input.genome} -m genome -l {wildcards.lineage} -o {wildcards.lineage} \\
              --out_path {params.outpath} -c {threads} -f {params.offline}
        cp {params.outpath}/{wildcards.lineage}/short_summary.specific.*.txt {output}
        """, env="qc", versions="busco --version")


rule fastp:
    input:
        lambda w: [x for x in (P.short_sample(w.sample)["r1"], P.short_sample(w.sample)["r2"]) if x],
    output:
        r1=O("00_Data/rna_short/clean/{sample}.R1.fq.gz"),
        r2=O("00_Data/rna_short/clean/{sample}.R2.fq.gz"),
        json=O("00_Data/qc/fastp/{sample}.json"),
        html=O("00_Data/qc/fastp/{sample}.html"),
    params:
        script=script("00_Data", "05_fastp.{sample}"),
        reads=lambda w: (f"-i {P.short_sample(w.sample)['r1']} -I {P.short_sample(w.sample)['r2']}"
                         if P.short_sample(w.sample)["r2"] else f"-i {P.short_sample(w.sample)['r1']}"),
        outs=lambda w, output: (f"-o {output.r1} -O {output.r2}" if P.short_sample(w.sample)["r2"]
                                else f"-o {output.r1}"),
        se_placeholder=lambda w, output: "" if P.short_sample(w.sample)["r2"] else f": | gzip > {output.r2}  # single-end: empty R2 placeholder",
    threads: 8
    log:
        **logs("00_Data", "05_fastp.{sample}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        fastp {params.reads} {params.outs} -w {threads} -j {output.json} -h {output.html}
        {params.se_placeholder}
        """, env="core", versions="fastp --version")


rule nanoplot:
    input:
        lambda w: P.long_sample(w.sample)["reads"],
    output:
        O("00_Data/qc/nanoplot/{sample}/NanoStats.txt"),
    params:
        script=script("00_Data", "06_nanoplot.{sample}"),
        outdir=O("00_Data/qc/nanoplot/{sample}"),
        fmt=lambda w: "--fasta" if any(P.long_sample(w.sample)["reads"].endswith(x)
                                       for x in (".fa", ".fasta", ".fa.gz", ".fasta.gz")) else "--fastq",
    threads: 4
    log:
        **logs("00_Data", "06_nanoplot.{sample}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        NanoPlot {params.fmt} {input} -o {params.outdir} -t {threads}
        """, env="core", versions="NanoPlot --version")


HOMOLOG_PC_ONLY = bool(config.get("homolog_protein_coding_only", True))


def homolog_gff(name):
    """Annotation used for the homolog downstream (gffread, GeMoMa)."""
    return O(f"00_Data/homologs/{name}/{name}.coding.gff3")


rule homolog_clean:
    """Keep only protein-coding genes of a homolog annotation (homolog_protein_coding_only: true);
    otherwise pass the annotation through unchanged."""
    input:
        lambda w: P.homolog(w.homolog)["gff3"],
    output:
        gff=O("00_Data/homologs/{homolog}/{homolog}.coding.gff3"),
        ids=O("00_Data/homologs/{homolog}/{homolog}.coding.ids"),
        report=O("00_Data/qc/homologs/{homolog}.clean_report.tsv"),
    params:
        script=script("00_Data", "07_homolog_clean.{homolog}"),
        mode="" if HOMOLOG_PC_ONLY else "--passthrough",
    log:
        **logs("00_Data", "07_homolog_clean.{homolog}"),
    shell:
        SH(PY("homolog_clean.py") + """ gff --gff {input} --out-gff {output.gff} --out-ids {output.ids} \\
            --out-report {output.report} {params.mode}
        cat {output.report}
        """)


rule homolog_prep:
    """Protein / CDS FASTA of the (cleaned) homolog annotation. A supplied protein FASTA is used
    when >= 50% of its records match kept transcripts (restricted to them); otherwise the proteins
    translated from the cleaned annotation are used."""
    input:
        genome=lambda w: P.homolog(w.homolog)["genome"],
        gff=rules.homolog_clean.output.gff,
        ids=rules.homolog_clean.output.ids,
    output:
        pep=O("00_Data/homologs/{homolog}/{homolog}.prot.faa"),
        cds=O("00_Data/homologs/{homolog}/{homolog}.cds.fa"),
    params:
        script=script("00_Data", "08_homolog_prep.{homolog}"),
        tmp=O("00_Data/homologs/{homolog}/{homolog}.gffread.faa"),
        link=O("00_Data/homologs/{homolog}/{homolog}.genome.fa"),
        given=lambda w: P.homolog(w.homolog)["protein"] or "",
        clean=PY("clean_pep.py"),
        select=PY("homolog_clean.py"),
    log:
        **logs("00_Data", "08_homolog_prep.{homolog}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        # gffread writes a .fai next to the FASTA: read the genome through a link in the results dir
        ln -sf {input.genome} {params.link}
        gffread -g {params.link} -x {output.cds} -y {params.tmp} {input.gff}
        {params.clean} {params.tmp} {output.pep}
        if [ -n "{params.given}" ]; then
            if {params.select} fasta --fasta {params.given} --ids {input.ids} --out {output.pep}.given; then
                mv {output.pep}.given {output.pep}
                echo "using supplied proteins restricted to protein-coding transcripts"
            else
                rm -f {output.pep}.given
                echo "WARNING: supplied protein IDs do not match the annotation; using gffread proteins" >&2
            fi
        fi
        rm -f {params.tmp}
        """, env="core", versions="gffread --version")


rule agat_config:
    """AGAT configuration used by every AGAT call. The progress bar is switched off: AGAT 1.7.0 aborts
    with "progress bar already finished" when its line count is exceeded (seen with GeMoMa GFF3).
    AGAT's own log files are off as well (they were written into the working directory)."""
    output:
        AGAT_CONFIG,
    params:
        script=script("00_Data", "00_agat_config"),
        tmp=O("00_Data/agat_config.tmp"),
    log:
        **logs("00_Data", "00_agat_config"),
    conda:
        "../envs/agat.yaml"
    shell:
        SH("""
        rm -rf {params.tmp} && mkdir -p {params.tmp} && cd {params.tmp}
        agat config --expose
        # log: false -> no agat_log_* folders in the working directory (AGAT output is in our logs/)
        sed -E 's/^progress_bar:.*/progress_bar: false/; s/^log:.*/log: false/' agat_config.yaml > {output}
        cd .. && rm -rf {params.tmp}
        grep -nE '^(progress_bar|log):' {output}
        """, env="agat")


rule homolog_stats:
    input:
        gff=lambda w: P.homolog(w.homolog)["gff3"],
        agat_cfg=AGAT_CONFIG,
    output:
        O("00_Data/qc/homologs/{homolog}.agat_stats.txt"),
    params:
        script=script("00_Data", "09_homolog_stats.{homolog}"),
    log:
        **logs("00_Data", "09_homolog_stats.{homolog}"),
    conda:
        "../envs/agat.yaml"
    shell:
        SH("""
        # AGAT names its temp folder after the input file in the current directory: one folder per job
        AGAT_WD={params.script}.agat; rm -rf "$AGAT_WD"; mkdir -p "$AGAT_WD"; cd "$AGAT_WD"
        agat_sp_statistics.pl --gff {input.gff} -o {output} --config {input.agat_cfg}
        cd - > /dev/null; rm -rf "$AGAT_WD"
        """, env="agat")


ORTHODB_SUBSET = O("00_Data/protein_db/orthodb_subset.faa")


rule orthodb_subset:
    """Proteins of the OrthoDB levels in orthodb_subset.levels (e.g. Brassicales), extracted from the
    OrthoDB partition on every run; keeps BRAKER4's protein copies small."""
    input:
        fasta=P.orthodb,
        levels=P.orthodb_levels_tab or [],
        l2s=P.orthodb_level2species_tab or [],
    output:
        fasta=ORTHODB_SUBSET,
        report=O("00_Data/qc/orthodb_subset.tsv"),
    params:
        script=script("00_Data", "09a_orthodb_subset"),
        select=" ".join(f"'{x}'" for x in P.orthodb_levels),
    log:
        **logs("00_Data", "09a_orthodb_subset"),
    shell:
        SH(PY("orthodb_subset.py") + """ --fasta {input.fasta} --levels {input.levels} \
            --level2species {input.l2s} --select {params.select} --out {output.fasta} --report {output.report}
        """)


rule protein_db:
    input:
        pep=expand(O("00_Data/homologs/{homolog}/{homolog}.prot.faa"), homolog=[h["name"] for h in P.homologs]),
        cds=expand(O("00_Data/homologs/{homolog}/{homolog}.cds.fa"), homolog=[h["name"] for h in P.homologs]),
        orthodb=[ORTHODB_SUBSET] if P.orthodb and P.orthodb_levels else [P.orthodb] if P.orthodb else [],
    output:
        all=PROTEINS,
        homolog=HOMOLOG_PROTEINS,
        cds=HOMOLOG_CDS,
    params:
        script=script("00_Data", "09b_protein_db"),
    log:
        **logs("00_Data", "09b_protein_db"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        # one input file -> symbolic link instead of a copy
        merge() {{ local out=$1; shift; if [ $# -eq 1 ]; then ln -sfn "$(readlink -f "$1")" "$out"; else cat "$@" > "$out"; fi; }}
        merge {output.homolog} {input.pep}
        merge {output.cds} {input.cds}
        orthodb="{input.orthodb}"
        if [ -n "$orthodb" ] && ! gzip -t "$orthodb" 2> /dev/null && [ ! -s {output.homolog} ]; then
            ln -sfn "$orthodb" {output.all}      # OrthoDB only, uncompressed: link
        else
            # OrthoDB (may be .gz) + homolog proteins, exact duplicates removed: a merged file is needed
            seqkit seq -w 60 $orthodb {output.homolog} | seqkit rmdup -s > {output.all}
        fi
        seqkit stats -T {output.all} {output.homolog}
        """, env="core")


rule data_qc:
    """Target that runs all raw-data assessments in 00_Data/qc."""
    input:
        O("00_Data/qc/seqkit/genome.stats.tsv"),
        O("00_Data/qc/assembly_quality.tsv"),
        expand(O("00_Data/qc/busco_genome/{lineage}/short_summary.txt"), lineage=P.busco_lineages),
        expand(O("00_Data/qc/fastp/{sample}.json"), sample=P.short_names),
        expand(O("00_Data/qc/nanoplot/{sample}/NanoStats.txt"), sample=P.long_names),
        expand(O("00_Data/qc/homologs/{homolog}.agat_stats.txt"), homolog=[h["name"] for h in P.homologs]),
        expand(O("00_Data/qc/homologs/{homolog}.clean_report.tsv"), homolog=[h["name"] for h in P.homologs]),
