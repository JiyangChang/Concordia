# 00_Data/qc: genome assembly quality assessment.
#   default : contiguity (scaffold + contig level), BUSCO genome mode (data.smk), LAI (after EDTA)
#   optional: Merqury (needs WGS reads), tidk telomeres (needs motif), NCBI FCS-GX (needs db + tax_id)
# Everything runs on the assembly as supplied (original sequence IDs, before renaming/filtering).

AQ = config.get("assembly_qc", {}) or {}
AQ_DIR = O("00_Data/qc/assembly")


rule assembly_stats:
    input:
        GENOME_RAW,
    output:
        summary=f"{AQ_DIR}/assembly_stats.tsv",
        seqs=f"{AQ_DIR}/sequences.tsv",
    params:
        script=script("00_Data", "10_assembly_stats"),
        gap=AQ.get("gap_min_n", 10),
    log:
        **logs("00_Data", "10_assembly_stats"),
    shell:
        SH(PY("assembly_stats.py") + """ --fasta {input} --gap-min-n {params.gap} \\
            --out-summary {output.summary} --out-seqs {output.seqs}
        cat {output.summary}
        """)


rule lai:
    """LTR Assembly Index from EDTA's intact LTRs + whole-genome RepeatMasker .out (EDTA library)."""
    input:
        genome=GENOME,
        edta=rules.edta.output.anno,          # EDTA finished -> intact LTR pass.list exists
        rmout=O("01_Repeat/RepeatMasker/genome.fa.out"),
    output:
        f"{AQ_DIR}/lai/genome.LAI",
    params:
        script=script("00_Data", "11_lai"),
        workdir=f"{AQ_DIR}/lai",
        edta_dir=O("01_Repeat/EDTA"),
        mono="-mono" if AQ.get("lai_mono", False) else "",
    threads: 8
    log:
        **logs("00_Data", "11_lai"),
    conda:
        "../envs/repeat.yaml"
    shell:
        SH("""
        # EDTA does not document the exact location of the intact-LTR pass.list; search known layouts.
        # Skip *.nmtf.pass.list (intact LTRs without the TG..CA motif): LAI uses the main list.
        PASS=$(ls {params.edta_dir}/*EDTA.raw/LTR/*.pass.list {params.edta_dir}/*EDTA.raw/*.pass.list 2>/dev/null \
               | grep -v '[.]nmtf[.]pass[.]list$' | head -n 1 || true)
        if [ -z "$PASS" ]; then echo "ERROR: no *.pass.list found under {params.edta_dir}" >&2; exit 1; fi
        echo "intact LTR list: $PASS"
        ln -sf {input.genome} genome.fa
        ln -sf "$PASS" genome.fa.pass.list
        ln -sf {input.rmout} genome.fa.out
        LAI -genome genome.fa -intact genome.fa.pass.list -all genome.fa.out -t {threads} {params.mono} 2>&1 | tee lai.log
        if [ -s genome.fa.out.LAI ]; then
            cp genome.fa.out.LAI {output}
        else
            # LAI exits 0 without a result when it is not applicable (e.g. LTR content < 5%)
            REASON=$(grep -E 'Error|not applicable' lai.log | sed -E 's/^[[:space:]]*//; s/【Error】 *//' | tr '\\n' ' ' || true)
            echo "# LAI not applicable: ${{REASON:-no result written by LAI}}" > {output}
            cat {output}
        fi
        """, env="repeat", workdir="{params.workdir}")


MQ = AQ.get("merqury", {}) or {}


rule merqury:
    input:
        genome=GENOME_RAW,
        reads=P.merqury_reads,
    output:
        qv=f"{AQ_DIR}/merqury/genome.qv",
        completeness=f"{AQ_DIR}/merqury/genome.completeness.stats",
    params:
        script=script("00_Data", "12_merqury"),
        workdir=f"{AQ_DIR}/merqury",
        k=MQ.get("k", 21),
        mem_gb=max(4, res("merqury")["mem_mb"] // 1000 - 4),
    threads: res("merqury")["threads"]
    resources:
        mem_mb=res("merqury")["mem_mb"],
        runtime=res("merqury")["runtime"],
    log:
        **logs("00_Data", "12_merqury"),
    conda:
        "../envs/merqury.yaml"
    shell:
        SH("""
        rm -rf reads.meryl
        meryl count k={params.k} threads={threads} memory={params.mem_gb} output reads.meryl {input.reads}
        ln -sf {input.genome} genome.fa
        merqury.sh reads.meryl genome.fa genome
        """, env="merqury", workdir="{params.workdir}", versions="meryl --version")


TEL = AQ.get("telomere", {}) or {}


rule tidk:
    input:
        genome=GENOME_RAW,
        seqs=rules.assembly_stats.output.seqs,
    output:
        windows=f"{AQ_DIR}/telomere/windows.tsv",
        table=f"{AQ_DIR}/telomere/telomere_per_sequence.tsv",
        summary=f"{AQ_DIR}/telomere/telomere_summary.tsv",
    params:
        script=script("00_Data", "13_tidk"),
        workdir=f"{AQ_DIR}/telomere",
        motif=P.telomere_motif or "",
        window=TEL.get("window", 10000),
        min_repeats=TEL.get("min_repeats", 50),
        summarise=PY("telomere_summary.py"),
    log:
        **logs("00_Data", "13_tidk"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        tidk search --string {params.motif} --window {params.window} --output tidk --dir . --extension tsv {input.genome}
        OUTTSV=$(ls tidk*.tsv 2>/dev/null | head -n 1 || true)
        if [ -z "$OUTTSV" ]; then echo "ERROR: tidk wrote no TSV" >&2; exit 1; fi
        cp "$OUTTSV" {output.windows}
        {params.summarise} --windows {output.windows} --seqs {input.seqs} --motif {params.motif} \\
            --min-repeats {params.min_repeats} --out {output.table} --out-summary {output.summary}
        cat {output.summary}
        """, env="core", workdir="{params.workdir}", versions="tidk --version")


FCS = AQ.get("fcs_gx", {}) or {}


rule fcs_gx:
    """NCBI FCS-GX contamination screen. Report only: no sequence is removed automatically."""
    input:
        genome=GENOME_RAW,
    output:
        f"{AQ_DIR}/fcs_gx/fcs_gx_report.txt",
    params:
        script=script("00_Data", "14_fcs_gx"),
        workdir=f"{AQ_DIR}/fcs_gx",
        cmd=FCS.get("cmd", "run_gx.py"),
        db=P.fcs_db or "",
        taxid=P.fcs_taxid or "",
    threads: res("fcs_gx")["threads"]
    resources:
        mem_mb=res("fcs_gx")["mem_mb"],
        runtime=res("fcs_gx")["runtime"],
    log:
        **logs("00_Data", "14_fcs_gx"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        ln -sf {input.genome} genome.fa
        export GX_NUM_CORES={threads}
        {params.cmd} --fasta genome.fa --out-dir . --gx-db {params.db} --tax-id {params.taxid}
        REP=$(ls *.fcs_gx_report.txt 2>/dev/null | head -n 1 || true)
        if [ -z "$REP" ]; then echo "ERROR: FCS-GX wrote no *.fcs_gx_report.txt" >&2; exit 1; fi
        cp "$REP" {output}
        """, env="core", workdir="{params.workdir}")


def assembly_qc_inputs():
    d = {"stats": rules.assembly_stats.output.summary,
         "busco": expand(O("00_Data/qc/busco_genome/{lineage}/short_summary.txt"), lineage=P.busco_lineages)}
    if P.run_lai:
        d["lai"] = rules.lai.output
    if P.run_merqury:
        d["qv"] = rules.merqury.output.qv
        d["completeness"] = rules.merqury.output.completeness
    if P.run_telomere:
        d["telomere"] = rules.tidk.output.summary
    if P.run_fcs_gx:
        d["fcs"] = rules.fcs_gx.output
    return d


def _opt(flag, key):
    return lambda w, input: f"{flag} {input[key]}" if key in input.keys() else ""


rule assembly_qc_summary:
    input:
        **assembly_qc_inputs(),
    output:
        O("00_Data/qc/assembly_quality.tsv"),
    params:
        script=script("00_Data", "15_assembly_qc_summary"),
        busco=" ".join(f"{lin}={O('00_Data/qc/busco_genome', lin, 'short_summary.txt')}" for lin in P.busco_lineages),
        lai=_opt("--lai", "lai"),
        qv=_opt("--merqury-qv", "qv"),
        comp=_opt("--merqury-completeness", "completeness"),
        tel=_opt("--telomere", "telomere"),
        fcs=_opt("--fcs", "fcs"),
    log:
        **logs("00_Data", "15_assembly_qc_summary"),
    shell:
        SH(PY("assembly_qc_summary.py") + """ --stats {input.stats} --busco {params.busco} \\
            {params.lai} {params.qv} {params.comp} {params.tel} {params.fcs} --out {output}
        cat {output}
        """)
