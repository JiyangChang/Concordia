# 02_RNA: short/long read alignment, merged BAMs, StringTie transcripts, EviAnn.


rule hisat2_index:
    input:
        GENOME,
    output:
        O("02_RNA/align_short/index/genome.1.ht2"),
    params:
        script=script("02_RNA", "01_hisat2_index"),
        prefix=O("02_RNA/align_short/index/genome"),
    threads: 16
    log:
        **logs("02_RNA", "01_hisat2_index"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        hisat2-build -p {threads} {input} {params.prefix}
        """, env="core", versions="hisat2 --version | head -n 1")


rule hisat2_align:
    input:
        idx=rules.hisat2_index.output,
        r1=O("00_Data/rna_short/clean/{sample}.R1.fq.gz"),
        r2=O("00_Data/rna_short/clean/{sample}.R2.fq.gz"),
    output:
        bam=O("02_RNA/align_short/{sample}.bam"),
        bai=O("02_RNA/align_short/{sample}.bam.bai"),
        summary=O("02_RNA/qc/{sample}.hisat2.txt"),
    params:
        script=script("02_RNA", "02_hisat2.{sample}"),
        prefix=O("02_RNA/align_short/index/genome"),
        reads=lambda w, input: (f"-1 {input.r1} -2 {input.r2}" if P.short_sample(w.sample)["r2"] else f"-U {input.r1}"),
        strand=lambda w: HISAT2_STRAND[P.short_sample(w.sample)["strandness"]],
        max_intron=P.max_intron,
    threads: 16
    resources:
        mem_mb=32000,
    log:
        **logs("02_RNA", "02_hisat2.{sample}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        hisat2 -p {threads} --dta -x {params.prefix} {params.reads} {params.strand} \\
               --max-intronlen {params.max_intron} --summary-file {output.summary} \\
          | samtools sort -@ 4 -o {output.bam} -
        samtools index {output.bam}
        """, env="core")


rule minimap2_align:
    input:
        genome=GENOME,
        reads=lambda w: P.long_sample(w.sample)["reads"],
    output:
        bam=O("02_RNA/align_long/{sample}.bam"),
        bai=O("02_RNA/align_long/{sample}.bam.bai"),
    params:
        script=script("02_RNA", "03_minimap2.{sample}"),
        preset=lambda w: MINIMAP2_PRESET[P.long_sample(w.sample)["platform"]],
        max_intron=P.max_intron,
    threads: 16
    resources:
        mem_mb=32000,
    log:
        **logs("02_RNA", "03_minimap2.{sample}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        minimap2 -t {threads} {params.preset} -G {params.max_intron} {input.genome} {input.reads} \\
          | samtools sort -@ 4 -o {output.bam} -
        samtools index {output.bam}
        """, env="core", versions="minimap2 --version")


rule merge_bam:
    input:
        lambda w: P.merged_bam_members(w.group),
    output:
        bam=O("02_RNA/merged/{group}.bam"),
        bai=O("02_RNA/merged/{group}.bam.bai"),
        flagstat=O("02_RNA/qc/{group}.flagstat.txt"),
    params:
        script=script("02_RNA", "04_merge_bam.{group}"),
    threads: 8
    log:
        **logs("02_RNA", "04_merge_bam.{group}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        samtools merge -f -@ {threads} -o {output.bam} {input}
        samtools index {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat}
        """, env="core", versions="samtools --version | head -n 1")


rule stringtie:
    input:
        O(f"02_RNA/merged/{TX_GROUP}.bam"),
    output:
        O(f"02_RNA/stringtie/{TX_SOURCE}.gtf"),
    params:
        script=script("02_RNA", f"05_stringtie.{TX_SOURCE}"),
        long="-L" if TX_LONG else "",
    threads: 8
    log:
        **logs("02_RNA", f"05_stringtie.{TX_SOURCE}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        stringtie {params.long} -p {threads} -o {output} {input}
        """, env="core", versions="stringtie --version")


EVIANN_BAMS, EVIANN_TAG = P.eviann_reads()


rule eviann_reads:
    input:
        EVIANN_BAMS,
    output:
        O("02_RNA/EviAnn/reads.txt"),
    params:
        script=script("02_RNA", "06_eviann_reads"),
        line=" ".join(EVIANN_BAMS + [EVIANN_TAG]),
    log:
        **logs("02_RNA", "06_eviann_reads"),
    shell:
        SH("""
        # EviAnn tag decided from the data: {EVIANN_TAG} (only-ONT is provisionally bam_isoseq)
        printf '%s\\n' '{params.line}' > {output}
        """.replace("{EVIANN_TAG}", EVIANN_TAG))


EV = config.get("eviann", {}) or {}
_sp = EV.get("swissprot")
EVIANN_SWISSPROT = os.path.abspath(os.path.expanduser(str(_sp))) if _sp not in (None, "", "null") else None


rule eviann:
    """RNA-based annotation. EviAnn always needs a protein set to call CDS: with eviann.swissprot it gets
    UniProt-SwissProt (-s, no homolog data); otherwise it falls back to the homolog proteins (-p)."""
    input:
        genome=GENOME,
        reads=rules.eviann_reads.output,
        proteins=[EVIANN_SWISSPROT] if EVIANN_SWISSPROT else [HOMOLOG_PROTEINS],
        bams=EVIANN_BAMS,
    output:
        gff=O("02_RNA/EviAnn/genome.fa.pseudo_label.gff"),
    params:
        script=script("02_RNA", "07_eviann"),
        workdir=O("02_RNA/EviAnn"),
        max_intron=P.max_intron,
        prot_opt="-s" if EVIANN_SWISSPROT else "-p",
        prot_note=("UniProt-SwissProt (eviann.swissprot); homolog proteins are not used" if EVIANN_SWISSPROT
                   else "eviann.swissprot not set: homolog proteins (-p) are used"),
    threads: res("eviann")["threads"]
    resources:
        mem_mb=res("eviann")["mem_mb"],
        runtime=res("eviann")["runtime"],
    log:
        **logs("02_RNA", "07_eviann"),
    conda:
        "../envs/eviann.yaml"
    shell:
        SH("""
        ln -sf {input.genome} genome.fa
        # protein set: {params.prot_note}
        eviann.sh -t {threads} -g genome.fa -r {input.reads} {params.prot_opt} {input.proteins} -m {params.max_intron}
        """, env="eviann", workdir="{params.workdir}", versions="eviann.sh --version")


rule eviann_split:
    input:
        rules.eviann.output.gff,
    output:
        coding=O("02_RNA/EviAnn/eviann.coding.gff3"),
        lnc=O("02_RNA/EviAnn/eviann.lncRNA.gff3"),
        pseudo=O("02_RNA/EviAnn/eviann.pseudogene.gff3"),
    params:
        script=script("02_RNA", "08_eviann_split"),
    log:
        **logs("02_RNA", "08_eviann_split"),
    shell:
        SH(PY("split_coding.py") + """ --gff {input} --coding {output.coding} --noncoding {output.lnc} \\
            --pseudogenes {output.pseudo}
        """)



rule eviann_stats:
    """Gene-model statistics of the EviAnn annotation (coding and lncRNA), for 07_Report."""
    input:
        coding=rules.eviann_split.output.coding,
        lnc=rules.eviann_split.output.lnc,
    output:
        stats=O("02_RNA/EviAnn/stats/eviann_stats.tsv"),
        dist=O("02_RNA/EviAnn/stats/eviann_distributions.tsv"),
    params:
        script=script("02_RNA", "11_eviann_stats"),
    log:
        **logs("02_RNA", "11_eviann_stats"),
    shell:
        SH(PY("annotation_stats.py") + """ --set coding={input.coding} --set lncRNA={input.lnc} \
            --out {output.stats} --out-dist {output.dist}
        cat {output.stats}
        """)

# ---------- alignment statistics (for 02_RNA/qc and 07_Report) ----------

rule bam_sample_stats:
    """samtools stats + number of spliced primary alignments (CIGAR contains N) for one sample."""
    input:
        O("02_RNA/align_{kind}/{sample}.bam"),
    output:
        stats=O("02_RNA/qc/{kind}/{sample}.samtools_stats.txt"),
        spliced=O("02_RNA/qc/{kind}/{sample}.spliced.txt"),
    wildcard_constraints:
        kind="short|long",
    params:
        script=script("02_RNA", "09_bam_stats.{kind}.{sample}"),
    threads: 4
    log:
        **logs("02_RNA", "09_bam_stats.{kind}.{sample}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        samtools stats -@ {threads} {input} | grep '^SN' | cut -f 2- > {output.stats}
        samtools view -@ {threads} -F 0x904 {input} | awk '$6 ~ /N/ {{n++}} END {{print n+0}}' > {output.spliced}
        """, env="core")


def rna_sample_entries():
    """kind:sample:platform for every RNA sample."""
    e = [f"short:{s['sample']}:illumina" for s in P.short]
    e += [f"long:{s['sample']}:{s['platform']}" for s in P.long]
    return e


rule rna_align_summary:
    input:
        stats=[O(f"02_RNA/qc/short/{s}.samtools_stats.txt") for s in P.short_names]
              + [O(f"02_RNA/qc/long/{s}.samtools_stats.txt") for s in P.long_names],
        spliced=[O(f"02_RNA/qc/short/{s}.spliced.txt") for s in P.short_names]
                + [O(f"02_RNA/qc/long/{s}.spliced.txt") for s in P.long_names],
        hisat2=[O(f"02_RNA/qc/{s}.hisat2.txt") for s in P.short_names],
        merged=[O(f"02_RNA/qc/{g}.flagstat.txt") for g in P.merged_bam_groups()],
    output:
        samples=O("02_RNA/qc/alignment_summary.tsv"),
        merged=O("02_RNA/qc/merged_bam_summary.tsv"),
    params:
        script=script("02_RNA", "10_rna_align_summary"),
        qc=O("02_RNA/qc"),
        samples=" ".join(rna_sample_entries()),
        merged=" ".join(f"'{g}={';'.join(u)}'" for g, u in P.bam_usage().items()),
    log:
        **logs("02_RNA", "10_rna_align_summary"),
    shell:
        SH(PY("rna_align_stats.py") + """ --qc-dir {params.qc} --samples {params.samples} \\
            --merged {params.merged} --out-samples {output.samples} --out-merged {output.merged}
        cat {output.samples} {output.merged}
        """)


# ---------- intron evidence for 05_EVM/refine ----------

LR_WEIGHT = bool((config.get("refine", {}) or {}).get("lr_weight_by_coverage", True))


rule lr_sample_chains:
    """Per-read intron chains of one long-read sample. Each read is weighted by its number of
    full-length reads: from the sample's Iso-Seq count table (rna_long.tsv column `counts`) if given,
    else from full_length_coverage in its FASTA/FASTQ header (Iso-Seq HQ transcripts), else 1.
    Done per sample because HQ transcript names (transcript/0, ...) repeat between Iso-Seq runs."""
    input:
        bam=O("02_RNA/align_long/{sample}.bam"),
        reads=lambda w: P.long_sample(w.sample)["reads"],
        counts=lambda w: [P.long_sample(w.sample)["counts"]] if P.long_sample(w.sample)["counts"] else [],
    output:
        chains=O("02_RNA/introns/lr_chains/{sample}.tsv"),
        summary=O("02_RNA/introns/lr_chains/{sample}.weights.tsv"),
    params:
        script=script("02_RNA", "12_lr_chains.{sample}"),
        weights=lambda w, input: ("" if not LR_WEIGHT else
                                  f"--counts {input.counts}" if input.counts else f"--reads {input.reads}"),
        platform=lambda w: P.long_sample(w.sample)["platform"],
    threads: 4
    log:
        **logs("02_RNA", "12_lr_chains.{sample}"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        samtools view -@ {threads} -F 0x904 {input.bam} | """ + PY("extract_introns.py") + """ chains {params.weights} \\
            --sample {wildcards.sample} --platform {params.platform} --summary {output.summary} --out {output.chains}
        cat {output.summary}
        """, env="core")


rule lr_read_chains:
    """Long-read intron chains of all samples (refine: whole-chain support, R3/R4)."""
    input:
        chains=[O(f"02_RNA/introns/lr_chains/{s}.tsv") for s in P.long_names],
        summaries=[O(f"02_RNA/introns/lr_chains/{s}.weights.tsv") for s in P.long_names],
    output:
        chains=O("02_RNA/introns/lr_read_chains.tsv"),
        summary=O("02_RNA/qc/long_read_weights.tsv"),
    params:
        script=script("02_RNA", "12b_lr_read_chains"),
    log:
        **logs("02_RNA", "12b_lr_read_chains"),
    shell:
        SH("""
        awk 'FNR > 1 || NR == 1' {input.chains} > {output.chains}
        awk 'FNR > 1 || NR == 1' {input.summaries} > {output.summary}
        cat {output.summary}
        """)


rule sr_junctions:
    """Short-read splice junction counts (refine fallback when there are no long reads)."""
    input:
        O("02_RNA/merged/sr.bam"),
    output:
        O("02_RNA/introns/sr_junctions.tsv"),
    params:
        script=script("02_RNA", "13_sr_junctions"),
    threads: 4
    log:
        **logs("02_RNA", "13_sr_junctions"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        samtools view -@ {threads} -F 0x904 {input} | awk '$6 ~ /N/' | """ + PY("extract_introns.py") + """ junctions --out {output}
        echo "junctions: $(($(wc -l < {output}) - 1))"
        """, env="core")
