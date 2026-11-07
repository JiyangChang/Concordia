# BUSCO (protein mode, user-specified lineages) + OMArk + AGAT statistics for every annotation
# source (EviAnn, BRAKER4, GeMoMa of every homologous species, merged GeMoMa, EVM, refine, final).
# One set of rules is generated per source.

QC_NAMES = {"eviann": "EviAnn", "braker4": "BRAKER4", "gemoma": "GeMoMa_merged", "evm": "EVM",
            "refine": "refine", "final": "final"}


def qc_name(key):
    return QC_NAMES.get(key, key.replace("gemoma_", "GeMoMa_", 1))


for _src, _qcdir in QC_DIRS.items():
    _qc = O(_qcdir)
    _step = _qcdir.split("/")[0]

    rule:
        name: f"qc_longest_{_src}"
        input:
            gff=QC_GFF[_src],
            genome=GENOME,
            agat_cfg=AGAT_CONFIG,
        output:
            gff=f"{_qc}/longest.gff3",
            pep=f"{_qc}/longest.pep.fa",
        params:
            script=script(_step, f"qc_01_longest.{_src}"),
            raw=f"{_qc}/longest.raw.pep.fa",
            clean=PY("clean_pep.py"),
        log:
            **logs(_step, f"qc_01_longest.{_src}"),
        conda:
            "../envs/agat.yaml"
        shell:
            SH("""
            # AGAT names its temp folder after the input file in the current directory: one folder per job
            AGAT_WD={params.script}.agat; rm -rf "$AGAT_WD"; mkdir -p "$AGAT_WD"; cd "$AGAT_WD"
            agat_sp_keep_longest_isoform.pl --gff {input.gff} -o {output.gff} --config {input.agat_cfg}
            cd - > /dev/null; rm -rf "$AGAT_WD"
            gffread -g {input.genome} -y {params.raw} {output.gff}
            {params.clean} {params.raw} {output.pep}
            """, env="agat")

    rule:
        name: f"qc_busco_{_src}"
        input:
            f"{_qc}/longest.pep.fa",
        output:
            f"{_qc}/busco/{{lineage}}/short_summary.txt",
        params:
            script=script(_step, f"qc_02_busco.{_src}.{{lineage}}"),
            outpath=f"{_qc}/busco",
            offline=busco_offline(),
        threads: res("busco")["threads"]
        resources:
            mem_mb=res("busco")["mem_mb"],
            runtime=res("busco")["runtime"],
        log:
            **logs(_step, f"qc_02_busco.{_src}.{{lineage}}"),
        conda:
            "../envs/qc.yaml"
        shell:
            SH("""
            busco -i {input} -m protein -l {wildcards.lineage} -o {wildcards.lineage} \\
                  --out_path {params.outpath} -c {threads} -f {params.offline}
            cp {params.outpath}/{wildcards.lineage}/short_summary.specific.*.txt {output}
            """, env="qc", versions="busco --version")

    rule:
        name: f"qc_omark_{_src}"
        input:
            pep=f"{_qc}/longest.pep.fa",
            db=P.omark_db,
        output:
            f"{_qc}/omark/omark.sum",
        params:
            script=script(_step, f"qc_03_omark.{_src}"),
            outdir=f"{_qc}/omark",
        threads: 8
        resources:
            mem_mb=res("omark")["mem_mb"],      # OMAmer loads the whole database (LUCA.h5: ~12 GB RAM)
        log:
            **logs(_step, f"qc_03_omark.{_src}"),
        conda:
            "../envs/omark.yaml"
        shell:
            SH("""
            omamer search --db {input.db} --query {input.pep} --out {params.outdir}/proteins.omamer --nthreads {threads}
            omark -f {params.outdir}/proteins.omamer -d {input.db} -o {params.outdir}/result
            cp {params.outdir}/result/*.sum {output}
            """, env="omark", versions="omark --version")

    rule:
        name: f"qc_agat_{_src}"
        input:
            gff=QC_GFF[_src],
            agat_cfg=AGAT_CONFIG,
        output:
            f"{_qc}/agat/stats.txt",
        params:
            script=script(_step, f"qc_04_agat.{_src}"),
        log:
            **logs(_step, f"qc_04_agat.{_src}"),
        conda:
            "../envs/agat.yaml"
        shell:
            SH("""
            # AGAT names its temp folder after the input file in the current directory: one folder per job
            AGAT_WD={params.script}.agat; rm -rf "$AGAT_WD"; mkdir -p "$AGAT_WD"; cd "$AGAT_WD"
            agat_sp_statistics.pl --gff {input.gff} -o {output} --config {input.agat_cfg}
            cd - > /dev/null; rm -rf "$AGAT_WD"
            """, env="agat")


rule qc_summary:
    input:
        busco=[O(d, "busco", lin, "short_summary.txt") for d in QC_DIRS.values() for lin in P.busco_lineages],
        omark=[O(d, "omark/omark.sum") for d in QC_DIRS.values()],
        agat=[O(d, "agat/stats.txt") for d in QC_DIRS.values()],
    output:
        O("05_EVM/qc/summary.tsv"),
    params:
        script=script("05_EVM", "qc_05_summary"),
        lineages=" ".join(P.busco_lineages),
        entries=" ".join(f"{s}={O(d)}" for s, d in QC_DIRS.items()),
    log:
        **logs("05_EVM", "qc_05_summary"),
    shell:
        SH(PY("qc_summary.py") + """ --out {output} --lineages {params.lineages} -- {params.entries}
        cat {output}
        """)


rule reference_compare:
    """Optional (qc.reference_gff): CDS-level agreement of every annotation set with a reference."""
    input:
        ref=P.reference_gff or [],
        map=SEQID_MAP,
        sets=[QC_GFF[k] for k in QC_DIRS],
    output:
        O("05_EVM/qc/reference_comparison.tsv"),
    params:
        script=script("05_EVM", "qc_06_reference_compare"),
        sets=" ".join(f"--set {qc_name(k)}={QC_GFF[k]}" for k in QC_DIRS),
    log:
        **logs("05_EVM", "qc_06_reference_compare"),
    shell:
        SH(PY("reference_compare.py") + """ --reference {input.ref} --seqid-map {input.map} {params.sets} \\
            --out {output}
        """)
