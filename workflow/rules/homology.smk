# 04_Homology: GeMoMa and miniprot, run separately for every homologous species.
#   GeMoMa:   RNA intron/coverage evidence extracted once (ERE) -> one GeMoMaPipeline per species
#             -> per-species predictions combined with GAF + AnnotationFinalizer (what a
#                multi-reference GeMoMaPipeline does internally) -> 04_Homology/gemoma.gff3 (to EVM)
#   miniprot: genome index built once -> one alignment per species -> IDs prefixed with the species
#             and concatenated -> 04_Homology/miniprot/miniprot.gff3 (to EVM as PROTEIN evidence)

GM = config.get("gemoma", {}) or {}
GEMOMA_BAM = P.gemoma_ere_bam()
HOMOLOG_NAMES = [h["name"] for h in P.homologs]
GEMOMA_ERE_DIR = O("04_Homology/GeMoMa/ERE")
GEMOMA_STRANDED = P.strandness in ("FR", "RF")


def gemoma_cmd(res_class="gemoma"):
    mem_gb = max(4, res(res_class)["mem_mb"] // 1000 - 2)
    # the bioconda GeMoMa wrapper already inserts "CLI" (verified with gemoma 1.9); passing CLI again
    # makes GeMoMa print its tool list and exit 0 without running
    return str(GM.get("cmd", "GeMoMa -Xmx{mem_gb}g")).replace("{mem_gb}", str(mem_gb))


def gemoma_species_gff(name):
    return O(f"04_Homology/GeMoMa/{name}/gemoma.{name}.gff3")


def gemoma_unfiltered(name):
    return O(f"04_Homology/GeMoMa/{name}/unfiltered_predictions_from_species_0.gff")


GEMOMA_ERE_OUT = ([f"{GEMOMA_ERE_DIR}/coverage_forward.bedgraph", f"{GEMOMA_ERE_DIR}/coverage_reverse.bedgraph"]
                  if GEMOMA_STRANDED else [f"{GEMOMA_ERE_DIR}/coverage.bedgraph"])


def gemoma_rna_args():
    """GeMoMaPipeline arguments that use the extracted RNA-seq evidence (r=EXTRACTED)."""
    if not GEMOMA_BAM:
        return "r=NO"
    if GEMOMA_STRANDED:
        cov = (f"coverage=STRANDED coverage_forward={GEMOMA_ERE_OUT[0]} coverage_reverse={GEMOMA_ERE_OUT[1]}")
    else:
        cov = f"coverage=UNSTRANDED coverage_unstranded={GEMOMA_ERE_OUT[0]}"
    return f"r=EXTRACTED introns={GEMOMA_ERE_DIR}/introns.gff {cov}"


rule gemoma_ere:
    """Extract introns and coverage from the RNA-seq BAM once; shared by all per-species runs."""
    input:
        bam=[GEMOMA_BAM] if GEMOMA_BAM else [],
    output:
        introns=f"{GEMOMA_ERE_DIR}/introns.gff",
        cov=GEMOMA_ERE_OUT,
    params:
        script=script("04_Homology", "01_gemoma_ere"),
        workdir=GEMOMA_ERE_DIR,
        cmd=gemoma_cmd("gemoma_ere"),   # own resource class: ERE needs far less than GeMoMa
        strand=GEMOMA_STRAND[P.strandness],
    resources:
        mem_mb=res("gemoma_ere")["mem_mb"],
    log:
        **logs("04_Homology", "01_gemoma_ere"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        {params.cmd} ERE s={params.strand} m={input.bam} c=true outdir=.
        ls -la
        """, env="core", workdir="{params.workdir}")


rule gemoma:
    """GeMoMaPipeline with ONE reference species."""
    input:
        genome=GENOME,
        ref_genome=lambda w: P.homolog(w.homolog)["genome"],
        ref_gff=lambda w: homolog_gff(w.homolog),
        ere=rules.gemoma_ere.output if GEMOMA_BAM else [],
    output:
        final=O("04_Homology/GeMoMa/{homolog}/final_annotation.gff"),
        unfiltered=O("04_Homology/GeMoMa/{homolog}/unfiltered_predictions_from_species_0.gff"),
    params:
        script=script("04_Homology", "02_gemoma.{homolog}"),
        workdir=O("04_Homology/GeMoMa/{homolog}"),
        cmd=gemoma_cmd(),
        rna=gemoma_rna_args(),
        max_intron=P.max_intron,
    threads: res("gemoma")["threads"]
    resources:
        mem_mb=res("gemoma")["mem_mb"],
        runtime=res("gemoma")["runtime"],
    log:
        **logs("04_Homology", "02_gemoma.{homolog}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        {params.cmd} GeMoMaPipeline threads={threads} outdir=. \\
            GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true \\
            t={input.genome} s=own i={wildcards.homolog} a={input.ref_gff} g={input.ref_genome} \\
            {params.rna} GeMoMa.m={params.max_intron}
        """, env="core", workdir="{params.workdir}")


rule gemoma_species_gff:
    input:
        rules.gemoma.output.final,
    output:
        O("04_Homology/GeMoMa/{homolog}/gemoma.{homolog}.gff3"),
    params:
        script=script("04_Homology", "03_gemoma_species_gff.{homolog}"),
    log:
        **logs("04_Homology", "03_gemoma_species_gff.{homolog}"),
    shell:
        SH("""
        ln -sfn final_annotation.gff {output}
        """)


rule gemoma_merge:
    """Combine the per-species predictions like a multi-reference GeMoMaPipeline does:
    GAF (filter, cluster, pick representatives; evidence = number of species) + AnnotationFinalizer."""
    input:
        genome=GENOME,
        unfiltered=[gemoma_unfiltered(h) for h in HOMOLOG_NAMES],
    output:
        gff=O("04_Homology/gemoma.gff3"),
    params:
        script=script("04_Homology", "04_gemoma_merge"),
        workdir=O("04_Homology/GeMoMa/merged"),
        cmd=gemoma_cmd(),
        # no GAF prefix (p=): the per-species runs (i=<species>) already start every ID with "<species>_"
        gaf_in=" ".join(f"g={gemoma_unfiltered(h)}" for h in HOMOLOG_NAMES),
    resources:
        mem_mb=res("gemoma")["mem_mb"],
    log:
        **logs("04_Homology", "04_gemoma_merge"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        {params.cmd} GAF {params.gaf_in} outdir=.
        {params.cmd} AnnotationFinalizer g={input.genome} a=filtered_predictions.gff rename=NO outdir=.
        cp final_annotation.gff {output.gff}
        """, env="core", workdir="{params.workdir}")


rule miniprot_index:
    input:
        genome=GENOME,
    output:
        O("04_Homology/miniprot/genome.mpi"),
    params:
        script=script("04_Homology", "05_miniprot_index"),
    threads: 16
    log:
        **logs("04_Homology", "05_miniprot_index"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        miniprot -t {threads} -d {output} {input.genome}
        """, env="core", versions="miniprot --version")


rule miniprot:
    """Proteins of ONE homologous species aligned to the genome; IDs prefixed with the species name."""
    input:
        index=rules.miniprot_index.output,
        proteins=O("00_Data/homologs/{homolog}/{homolog}.prot.faa"),
    output:
        gff=O("04_Homology/miniprot/species/{homolog}.gff3"),
    params:
        script=script("04_Homology", "06_miniprot.{homolog}"),
        max_intron=P.max_intron,
    threads: 16
    log:
        **logs("04_Homology", "06_miniprot.{homolog}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        # miniprot numbers its alignments MP000001... in every run: prefix them with the species
        miniprot -t {threads} --gff -G {params.max_intron} {input.index} {input.proteins} \\
          | sed -E 's/(ID|Parent)=MP/\\1={wildcards.homolog}_MP/g' > {output.gff}
        """, env="core", versions="miniprot --version")


rule miniprot_merge:
    input:
        expand(O("04_Homology/miniprot/species/{homolog}.gff3"), homolog=HOMOLOG_NAMES),
    output:
        O("04_Homology/miniprot/miniprot.gff3"),
    params:
        script=script("04_Homology", "07_miniprot_merge"),
    log:
        **logs("04_Homology", "07_miniprot_merge"),
    shell:
        SH("""
        {{ echo '##gff-version 3'; grep -hv '^##gff-version' {input} || true; }} > {output}
        """)


rule homology_stats:
    """Gene-model statistics of every per-species GeMoMa set and the merged set, per-species
    contribution, GeMoMa evidence attributes and per-species miniprot alignments, for 07_Report."""
    input:
        merged=rules.gemoma_merge.output.gff,
        species=[gemoma_species_gff(h) for h in HOMOLOG_NAMES],
        unfiltered=[gemoma_unfiltered(h) for h in HOMOLOG_NAMES],
        miniprot=expand(O("04_Homology/miniprot/species/{homolog}.gff3"), homolog=HOMOLOG_NAMES),
        proteins=expand(O("00_Data/homologs/{homolog}/{homolog}.prot.faa"), homolog=HOMOLOG_NAMES),
    output:
        stats=O("04_Homology/stats/gemoma_stats.tsv"),
        dist=O("04_Homology/stats/gemoma_distributions.tsv"),
        refs=O("04_Homology/stats/gemoma_references.tsv"),
        attrs=O("04_Homology/stats/gemoma_attributes.tsv"),
        miniprot=O("04_Homology/stats/miniprot_species.tsv"),
    params:
        script=script("04_Homology", "08_homology_stats"),
        sets=" ".join(f"--set GeMoMa_{h}={gemoma_species_gff(h)}" for h in HOMOLOG_NAMES),
        refs=" ".join(f"--reference {h}={gemoma_unfiltered(h)},{gemoma_species_gff(h)}" for h in HOMOLOG_NAMES),
        mp=" ".join(f"--species {h}={O('00_Data/homologs', h, h + '.prot.faa')},{O('04_Homology/miniprot/species', h + '.gff3')}"
                    for h in HOMOLOG_NAMES),
    log:
        **logs("04_Homology", "08_homology_stats"),
    shell:
        SH(PY("annotation_stats.py") + """ {params.sets} --set GeMoMa_merged={input.merged} \\
            --out {output.stats} --out-dist {output.dist}
        """ + PY("gemoma_evidence.py") + """ --final {input.merged} {params.refs} \\
            --out-references {output.refs} --out-attributes {output.attrs}
        """ + PY("miniprot_stats.py") + """ {params.mp} --out {output.miniprot}
        cat {output.stats} {output.refs} {output.attrs} {output.miniprot}
        """)
