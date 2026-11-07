# 05_EVM: EVM inputs -> EVM (standard v2.1.0 by default; Rust only if it supports --repeats)
#         -> refine against EviAnn -> filter -> rename -> restore sequence IDs -> export.

EVM = config.get("evm", {}) or {}
RF = config.get("refine", {}) or {}
SEG, OVL = P.evm_overlap()
EVM_SOURCE_GFF = {k: QC_GFF[k] for k, _ in GENEPRED}          # no EviAnn without RNA
EVM_SOURCE_NAME = dict(GENEPRED)


rule no_rna_placeholders:
    """Without RNA-seq (rna.allow_none): empty EviAnn coding / lncRNA GFF3, so that refine (every locus
    becomes R5), naming and export run unchanged."""
    output:
        coding=f"{NO_RNA_DIR}/eviann.coding.gff3",
        lnc=f"{NO_RNA_DIR}/eviann.lncRNA.gff3",
    shell:
        "printf '##gff-version 3\\n' > {output.coding}; printf '##gff-version 3\\n' > {output.lnc}"


rule agat_normalize:
    input:
        gff=lambda w: EVM_SOURCE_GFF[w.src],
        agat_cfg=AGAT_CONFIG,
    output:
        O("05_EVM/inputs/{src}.agat.gff3"),
    params:
        script=script("05_EVM", "01_agat_normalize.{src}"),
    log:
        **logs("05_EVM", "01_agat_normalize.{src}"),
    conda:
        "../envs/agat.yaml"
    shell:
        SH("""
        # AGAT names its temp folder after the input file in the current directory: one folder per job
        AGAT_WD={params.script}.agat; rm -rf "$AGAT_WD"; mkdir -p "$AGAT_WD"; cd "$AGAT_WD"
        agat_convert_sp_gxf2gxf.pl --gff {input.gff} -o {output} --config {input.agat_cfg}
        cd - > /dev/null; rm -rf "$AGAT_WD"
        """, env="agat", versions="agat_convert_sp_gxf2gxf.pl --help | head -n 5")


rule evm_inputs:
    input:
        genepred=expand(O("05_EVM/inputs/{src}.agat.gff3"), src=list(EVM_SOURCE_GFF)),
        miniprot=O("04_Homology/miniprot/miniprot.gff3"),
        transcripts=[TX_GTF] if HAS_RNA else [],
        repeats=O("01_Repeat/RepeatMasker/genome.fa.out.gff"),
    output:
        genepred=O("05_EVM/inputs/gene_predictions.gff3"),
        protein=O("05_EVM/inputs/protein_alignments.gff3"),
        transcript=O("05_EVM/inputs/transcript_alignments.gff3"),
        repeats=O("05_EVM/inputs/repeats.gff3"),
        weights=O("05_EVM/inputs/weights.txt"),
    params:
        script=script("05_EVM", "02_evm_inputs"),
        prep=PY("evm_prepare.py"),
        tmpdir=O("05_EVM/inputs/tmp"),
        weights=" ".join(l.replace("\t", ":") for l in P.evm_weights()),
        tx_source=TX_SOURCE or "-",
        genepred=" ".join(f"{k}:{n}" for k, n in GENEPRED),
        inputs=O("05_EVM/inputs"),
    log:
        **logs("05_EVM", "02_evm_inputs"),
    shell:
        SH("""
        mkdir -p {params.tmpdir}
        : > {output.genepred}
        for KN in {params.genepred}; do          # key:EVM source name, e.g. braker4:BRAKER4
            {params.prep} genepred --gff {params.inputs}/${{KN%%:*}}.agat.gff3 --source ${{KN#*:}} --out {params.tmpdir}/${{KN%%:*}}.gff3
            cat {params.tmpdir}/${{KN%%:*}}.gff3 >> {output.genepred}
        done
        {params.prep} miniprot --gff {input.miniprot} --out {output.protein}
        if [ -n "{input.transcripts}" ]; then
            {params.prep} transcripts --gtf {input.transcripts} --source {params.tx_source} --out {output.transcript}
        else
            : > {output.transcript}              # no RNA-seq: no TRANSCRIPT evidence
        fi
        {params.prep} repeats --gff {input.repeats} --out {output.repeats}
        {params.prep} weights --out {output.weights} {params.weights}
        cat {output.weights}
        """)


rule evm_impl_check:
    """Decide which EVM implementation to use. Rust is used only if its --help lists --repeats."""
    output:
        O("05_EVM/logs/evm_impl_check.txt"),
    params:
        script=script("05_EVM", "03_evm_impl_check"),
        want=EVM.get("implementation", "standard"),
        rust_bin=os.path.join(project_path(EVM.get("rust_repo", "external/EVidenceModeler_rust")),
                              "target/release/EVidenceModeler"),
    log:
        **logs("05_EVM", "03_evm_impl_check"),
    conda:
        "../envs/evm.yaml"
    shell:
        SH("""
        STD_BIN=$(command -v EVidenceModeler || true)
        if [ "{params.want}" = "rust" ]; then
            if [ -x "{params.rust_bin}" ] && "{params.rust_bin}" --help 2>&1 | grep -q -- '--repeats'; then
                printf 'rust\\t%s\\tRust build lists --repeats\\n' "{params.rust_bin}" > {output}
            else
                printf 'standard\\t%s\\tfallback: Rust binary missing or --repeats not listed in --help\\n' "$STD_BIN" > {output}
            fi
        else
            printf 'standard\\t%s\\tconfig evm.implementation=standard\\n' "$STD_BIN" > {output}
        fi
        if [ -z "$(cut -f2 {output})" ]; then echo "ERROR: EVidenceModeler not found in PATH" >&2; exit 1; fi
        cat {output}
        """, env="evm")


# ---------- EVM weight tuning (grid search on one sequence) ----------

TUNE = P.evm_tune
GRID = P.evm_grid() if TUNE["enabled"] else []
GRID_IDS = [c for c, _ in GRID]
TUNE_DIR = O("05_EVM/tune")
TUNE_INPUTS = ["gene_predictions.gff3", "protein_alignments.gff3", "transcript_alignments.gff3", "repeats.gff3"]


rule evm_tune_inputs:
    """EVM inputs restricted to the tuning region (evm.tune.region; auto = longest sequence)."""
    input:
        genome=GENOME,
        map=SEQID_MAP,
        gffs=[O("05_EVM/inputs", f) for f in TUNE_INPUTS],
    output:
        genome=f"{TUNE_DIR}/inputs/genome.fa",
        region=f"{TUNE_DIR}/inputs/region.tsv",
        gffs=[f"{TUNE_DIR}/inputs/{f}" for f in TUNE_INPUTS],
    params:
        script=script("05_EVM", "03a_evm_tune_inputs"),
        region=TUNE["region"],
        outdir=f"{TUNE_DIR}/inputs",
    log:
        **logs("05_EVM", "03a_evm_tune_inputs"),
    shell:
        SH(PY("subset_evm_inputs.py") + """ --genome {input.genome} --seqid-map {input.map} \\
            --region {params.region} --outdir {params.outdir} {input.gffs}
        cat {output.region}
        """)


rule evm_tune_weights:
    """One weights file per grid combination (c00 = configured evm.weights)."""
    output:
        [f"{TUNE_DIR}/grid/{c}/weights.txt" for c in GRID_IDS],
    params:
        script=script("05_EVM", "03b_evm_tune_weights"),
        cmds="\n".join(
            f"printf '%s\\n' " + " ".join(f"'{l}'" for l in P.evm_weight_lines(w)) + f" > {TUNE_DIR}/grid/{c}/weights.txt"
            for c, w in GRID),
    log:
        **logs("05_EVM", "03b_evm_tune_weights"),
    shell:
        SH("""
        {params.cmds}
        echo "grid combinations: """ + str(len(GRID_IDS)) + """"
        """)


rule evm_tune_run:
    """EVM on the tuning region with the weights of one grid combination."""
    input:
        check=rules.evm_impl_check.output,
        genome=rules.evm_tune_inputs.output.genome,
        gffs=rules.evm_tune_inputs.output.gffs,
        weights=f"{TUNE_DIR}/grid/{{combo}}/weights.txt",
    output:
        f"{TUNE_DIR}/grid/{{combo}}/tune.EVM.gff3",
    params:
        script=script("05_EVM", "03c_evm_tune_run.{combo}"),
        workdir=f"{TUNE_DIR}/grid/{{combo}}",
        d=f"{TUNE_DIR}/inputs",
        seg=SEG,
        ovl=OVL,
    threads: 4
    resources:
        mem_mb=8000,
    log:
        **logs("05_EVM", "03c_evm_tune_run.{combo}"),
    conda:
        "../envs/evm.yaml"
    shell:
        SH("""
        BIN=$(cut -f2 {input.check})
        CPU_ARG=""
        if "$BIN" --help 2>&1 | grep -q -- '--CPU'; then CPU_ARG="--CPU {threads}"; fi
        TX_ARG=""                                  # empty without RNA-seq
        if [ -s {params.d}/transcript_alignments.gff3 ]; then TX_ARG="--transcript_alignments {params.d}/transcript_alignments.gff3"; fi
        "$BIN" --sample_id tune --genome {input.genome} --weights {input.weights} \\
            --gene_predictions {params.d}/gene_predictions.gff3 \\
            --protein_alignments {params.d}/protein_alignments.gff3 \\
            $TX_ARG \\
            --repeats {params.d}/repeats.gff3 \\
            --segmentSize {params.seg} --overlapSize {params.ovl} $CPU_ARG
        [ -f tune.EVM.gff3 ] || : > tune.EVM.gff3
        """, env="evm", workdir="{params.workdir}")


rule evm_tune_score:
    """Score every combination against long reads / junctions and miniprot; pick the best."""
    input:
        unpack(lambda w: refine_support_input()),
        gffs=[f"{TUNE_DIR}/grid/{c}/tune.EVM.gff3" for c in GRID_IDS],
        weights=[f"{TUNE_DIR}/grid/{c}/weights.txt" for c in GRID_IDS],
        region=rules.evm_tune_inputs.output.region,
        miniprot=O("04_Homology/miniprot/miniprot.gff3"),
        stringtie=[TX_GTF] if HAS_RNA else [],
    output:
        table=f"{TUNE_DIR}/tuning_results.tsv",
        best=f"{TUNE_DIR}/best_weights.txt",
    params:
        script=script("05_EVM", "03d_evm_tune_score"),
        combos=" ".join(f"--combo {c}={TUNE_DIR}/grid/{c}/tune.EVM.gff3" for c in GRID_IDS),
        reads=lambda w, input: (f"--lr-chains {input.lr}" if "lr" in input.keys() else
                                f"--sr-junctions {input.sr}" if "sr" in input.keys() else ""),
        min_lr=RF.get("min_lr_support", 2),
        min_sr=RF.get("min_sr_junction_reads", 3),
        sr_ratio=RF.get("sr_min_junction_ratio", 0.1),
        sr_cov=RF.get("sr_stringtie_min_cov", 3),
    log:
        **logs("05_EVM", "03d_evm_tune_score"),
    shell:
        SH(PY("evm_tune_score.py") + """ {params.combos} --weights-dir """ + f"{TUNE_DIR}/grid" + """ \\
            --region-file {input.region} --miniprot {input.miniprot} {params.reads} \\
            --min-lr-support {params.min_lr} --min-sr-junction {params.min_sr} \\
            --stringtie "{input.stringtie}" --sr-min-junction-ratio {params.sr_ratio} --sr-stringtie-min-cov {params.sr_cov} \\
            --out {output.table} --best {output.best}
        head -n 6 {output.table}
        """)


EVM_WEIGHTS = (f"{TUNE_DIR}/best_weights.txt" if TUNE["enabled"] and TUNE["use_best"]
               else O("05_EVM/inputs/weights.txt"))


rule evm:
    input:
        genome=GENOME,
        check=rules.evm_impl_check.output,
        genepred=rules.evm_inputs.output.genepred,
        protein=rules.evm_inputs.output.protein,
        transcript=rules.evm_inputs.output.transcript,
        repeats=rules.evm_inputs.output.repeats,
        weights=EVM_WEIGHTS,
    output:
        gff=O("05_EVM/evm", f"{PREFIX}.EVM.gff3"),
    params:
        script=script("05_EVM", "04_evm"),
        workdir=O("05_EVM/evm"),
        seg=SEG,
        ovl=OVL,
    threads: res("evm")["threads"]
    resources:
        mem_mb=res("evm")["mem_mb"],
        runtime=res("evm")["runtime"],
    log:
        **logs("05_EVM", "04_evm"),
    conda:
        "../envs/evm.yaml"
    shell:
        SH("""
        BIN=$(cut -f2 {input.check})
        echo "EVM implementation: $(cut -f1,3 {input.check})"
        CPU_ARG=""
        if "$BIN" --help 2>&1 | grep -q -- '--CPU'; then CPU_ARG="--CPU {threads}"; fi
        TX_ARG=""                                  # empty without RNA-seq
        if [ -s {input.transcript} ]; then TX_ARG="--transcript_alignments {input.transcript}"; fi
        "$BIN" --sample_id """ + PREFIX + """ --genome {input.genome} --weights {input.weights} \\
            --gene_predictions {input.genepred} \\
            --protein_alignments {input.protein} \\
            $TX_ARG \\
            --repeats {input.repeats} \\
            --segmentSize {params.seg} --overlapSize {params.ovl} $CPU_ARG
        """, env="evm", workdir="{params.workdir}")


def refine_support_input():
    """Long-read intron chains if long reads exist, else short-read junctions (if enabled)."""
    if P.has_long:
        return {"lr": O("02_RNA/introns/lr_read_chains.tsv")}
    if RF.get("sr_junction_fallback", True) and P.has_short:
        return {"sr": O("02_RNA/introns/sr_junctions.tsv")}
    return {}


rule refine:
    """EVM vs EviAnn refinement (rules R1-R6, see workflow/scripts/refine.py)."""
    input:
        unpack(lambda w: refine_support_input()),
        evm=rules.evm.output.gff,
        eviann=EVIANN_CODING,                 # empty placeholder without RNA -> every locus is R5
        genome=GENOME,
        stringtie=[TX_GTF] if HAS_RNA else [],
        miniprot=O("04_Homology/miniprot/miniprot.gff3"),
        gemoma=QC_GFF["gemoma"],
        braker4=QC_GFF["braker4"],
        support=O("03_Ab_initio/braker4.gene_support.tsv"),
        te=O("01_Repeat/EDTA/genome.fa.mod.EDTA.TEanno.gff3"),
    output:
        gff=QC_GFF["refine"],
        decisions=O("05_EVM/refine/refine_decisions.tsv"),
    params:
        script=script("05_EVM", "05_refine"),
        reads=lambda w, input: (f"--lr-chains {input.lr}" if "lr" in input.keys() else
                                f"--sr-junctions {input.sr}" if "sr" in input.keys() else ""),
        min_aa=RF.get("min_cds_aa", 50),
        max_te=RF.get("max_te_overlap", 0.5),
        min_lr=RF.get("min_lr_support", 2),
        min_sr=RF.get("min_sr_junction_reads", 3),
        se_ov=RF.get("single_exon_overlap", 0.5),
        sr_ratio=RF.get("sr_min_junction_ratio", 0.1),
        sr_cov=RF.get("sr_stringtie_min_cov", 3),
        imp="true" if RF.get("import_utr_isoforms", True) else "false",
        rm="true" if RF.get("remove_unsupported_evm", False) else "false",
    resources:
        mem_mb=res("refine")["mem_mb"],
    log:
        **logs("05_EVM", "05_refine"),
    shell:
        SH(PY("refine.py") + """ --evm {input.evm} --eviann {input.eviann} --genome {input.genome} {params.reads} \
            --stringtie "{input.stringtie}" --proteins {input.miniprot} {input.gemoma} \
            --braker4-gff {input.braker4} --braker4-support {input.support} --te-gff {input.te} \
            --min-cds-aa {params.min_aa} --max-te-overlap {params.max_te} --min-lr-support {params.min_lr} \
            --min-sr-junction {params.min_sr} --single-exon-overlap {params.se_ov} \
            --sr-min-junction-ratio {params.sr_ratio} --sr-stringtie-min-cov {params.sr_cov} \
            --import-utr-isoforms {params.imp} --remove-unsupported-evm {params.rm} \
            --out {output.gff} --decisions {output.decisions}
        grep -v '^#' {output.decisions} | cut -f 8,9 | sort | uniq -c | sort -k1,1nr
        """)


rule final_pep:
    input:
        genome=GENOME,
        gff=QC_GFF["refine"],
    output:
        O("05_EVM/final/work/refined.pep.fa"),
    params:
        script=script("05_EVM", "06_final_pep"),
        raw=O("05_EVM/final/work/refined.raw.pep.fa"),
    log:
        **logs("05_EVM", "06_final_pep"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        gffread -g {input.genome} -y {params.raw} {input.gff}
        """ + PY("clean_pep.py") + """ {params.raw} {output}
        """, env="core")


rule final_tesorter:
    input:
        rules.final_pep.output,
    output:
        O("05_EVM/final/work/tesorter/refined.cls.tsv"),
    params:
        script=script("05_EVM", "07_final_tesorter"),
        workdir=O("05_EVM/final/work/tesorter"),
    threads: 8
    log:
        **logs("05_EVM", "07_final_tesorter"),
    conda:
        "../envs/repeat.yaml"
    shell:
        SH("""
        TEsorter {input} -st prot -p {threads} -pre refined
        [ -f refined.cls.tsv ] || : > refined.cls.tsv
        """, env="repeat", workdir="{params.workdir}")


rule final_filter:
    """Remove a gene for TE reasons only if CDS-TE overlap > max_te_overlap AND it has no RNA
    evidence (StringTie / EviAnn / BRAKER4 gene_support); also remove CDS < min_cds_aa.
    TEsorter hits are only recorded in filter_decisions.tsv."""
    input:
        gff=QC_GFF["refine"],
        te=O("01_Repeat/EDTA/genome.fa.mod.EDTA.TEanno.gff3"),
        transcripts=[TX_GTF] if HAS_RNA else [],
        decisions=O("05_EVM/refine/refine_decisions.tsv"),
        miniprot=O("04_Homology/miniprot/miniprot.gff3"),
        gemoma=QC_GFF["gemoma"],
        braker4=QC_GFF["braker4"],
        support=O("03_Ab_initio/braker4.gene_support.tsv"),
        tesorter=rules.final_tesorter.output,
    output:
        gff=O("05_EVM/final/work/filtered.gff3"),
        removed=O("05_EVM/final/removed_genes.tsv"),
        decisions=O("05_EVM/final/filter_decisions.tsv"),
    params:
        script=script("05_EVM", "08_final_filter"),
        min_aa=RF.get("min_cds_aa", 50),
        max_te=RF.get("max_te_overlap", 0.5),
        # without RNA-seq, protein evidence (miniprot / GeMoMa) replaces RNA evidence for the TE criterion
        prot=lambda w, input: "" if HAS_RNA else f"--protein-evidence {input.miniprot} {input.gemoma}",
    log:
        **logs("05_EVM", "08_final_filter"),
    shell:
        SH(PY("filter_final.py") + """ --gff {input.gff} \\
            --te-gff {input.te} --max-te-overlap {params.max_te} \\
            --transcripts "{input.transcripts}" --refine-decisions {input.decisions} {params.prot} \\
            --braker4-gff {input.braker4} --braker4-support {input.support} \\
            --tesorter {input.tesorter} --min-cds-aa {params.min_aa} \\
            --out {output.gff} --removed {output.removed} --decisions {output.decisions}
        """)


NM = config.get("naming", {}) or {}
LNC_INTERNAL = O("05_EVM/final/work/lncRNA.named.gff3")


rule final_names:
    """Chromosome-based gene names (e.g. Mol01G000010 / Mol01L000010 for lncRNA), UTR features,
    recomputed CDS phases, old -> new ID table. See workflow/scripts/rename_genes.py."""
    input:
        coding=rules.final_filter.output.gff,
        lnc=EVIANN_LNC,
        map=SEQID_MAP,
    output:
        coding=FINAL_INTERNAL_ALL,
        lnc=LNC_INTERNAL,
        ids=O("05_EVM/final/id_map.tsv"),
    params:
        script=script("05_EVM", "09_final_names"),
        scheme=NM.get("scheme", "chromosome"),
        step=NM.get("step", 10),
        digits=NM.get("digits", 6),
        regex=NM.get("chrom_regex", r"^(?:chr(?:omosome)?)?[_-]?0*(\d+)$"),
        unplaced=NM.get("unplaced_code", "00"),
        codes=" ".join(f"--code '{k}={v}'" for k, v in (NM.get("chrom_codes") or {}).items()),
        gtag=NM.get("coding_tag", "G"),
        ltag=NM.get("lncrna_tag", "L"),
    log:
        **logs("05_EVM", "09_final_names"),
    shell:
        SH(PY("rename_genes.py") + """ --coding {input.coding} --lncrna {input.lnc} --seqid-map {input.map} \\
            --prefix """ + PREFIX + """ --scheme {params.scheme} --step {params.step} --digits {params.digits} \\
            --chrom-regex '{params.regex}' --unplaced-code {params.unplaced} {params.codes} \\
            --coding-tag {params.gtag} --lncrna-tag {params.ltag} \\
            --out-coding {output.coding} --out-lncrna {output.lnc} --out-map {output.ids}
        """)


rule final_select_longest:
    """Main result = one transcript per gene: the longest CDS (ties: longest transcript, then the first).
    All isoforms are kept as a separate result (Mol.all_isoforms.*)."""
    input:
        FINAL_INTERNAL_ALL,
    output:
        gff=FINAL_INTERNAL,
        report=O("05_EVM/final/isoform_selection.tsv"),
    params:
        script=script("05_EVM", "09b_final_select_longest"),
    log:
        **logs("05_EVM", "09b_final_select_longest"),
    shell:
        SH(PY("select_longest.py") + """ --gff {input} --out {output.gff} --report {output.report}
        """)


rule final_export:
    """Restore the original sequence IDs and write GFF3 + CDS / protein / cDNA FASTA for the main
    result (one transcript per gene) and for the all-isoform result."""
    input:
        gff=FINAL_INTERNAL,
        all=FINAL_INTERNAL_ALL,
        lnc=LNC_INTERNAL,
        map=SEQID_MAP,
        genome=GENOME_RAW,
    output:
        gff=FINAL_GFF,
        cds=O("05_EVM/final", f"{PREFIX}.cds.fa"),
        pep=O("05_EVM/final", f"{PREFIX}.pep.fa"),
        cdna=O("05_EVM/final", f"{PREFIX}.cdna.fa"),
        all=FINAL_ALL_GFF,
        all_cds=O("05_EVM/final", f"{PREFIX}.all_isoforms.cds.fa"),
        all_pep=O("05_EVM/final", f"{PREFIX}.all_isoforms.pep.fa"),
        all_cdna=O("05_EVM/final", f"{PREFIX}.all_isoforms.cdna.fa"),
        lnc=O("05_EVM/final", f"{PREFIX}.lncRNA.gff3"),
        func_pep=FINAL_LONGEST_PEP,
        stops=O("05_EVM/final/internal_stops.tsv"),
    params:
        script=script("05_EVM", "10_final_export"),
        restore=PY("restore_seqids.py"),
        clean=PY("clean_pep.py"),
        mark=PY("mark_internal_stops.py"),
    log:
        **logs("05_EVM", "10_final_export"),
    conda:
        "../envs/core.yaml"
    shell:
        SH("""
        {params.restore} --gff {input.gff} --map {input.map} --out {output.gff}
        {params.restore} --gff {input.all} --map {input.map} --out {output.all}
        {params.restore} --gff {input.lnc} --map {input.map} --out {output.lnc}
        gffread -g {input.genome} -x {output.cds} -y {output.pep} -w {output.cdna} {output.gff}
        gffread -g {input.genome} -x {output.all_cds} -y {output.all_pep} -w {output.all_cdna} {output.all}
        # transcripts whose CDS translation has internal stop codons are kept and marked internal_stop=yes
        {params.mark} --gff {output.gff} --pep {output.pep} --report {output.stops}
        {params.mark} --gff {output.all} --pep {output.all_pep}
        {params.clean} {output.pep} {output.func_pep}     # proteins for InterProScan / eggNOG (internal stops -> X)
        """, env="core")


# ---------- statistics for 07_Report ----------

SOURCE_ARGS = " ".join(f"--source {n}={QC_GFF[k]}" for k, n in GENEPRED)


rule evm_stats:
    """EVM gene-model statistics + agreement of EVM transcripts with BRAKER4 / GeMoMa / EviAnn."""
    input:
        evm=QC_GFF["evm"],
        sources=[QC_GFF[k] for k, _ in GENEPRED],
    output:
        stats=O("05_EVM/stats/evm_stats.tsv"),
        dist=O("05_EVM/stats/evm_distributions.tsv"),
        support=O("05_EVM/stats/evm_source_support.tsv"),
        combos=O("05_EVM/stats/evm_source_combinations.tsv"),
    params:
        script=script("05_EVM", "12_evm_stats"),
        sources=SOURCE_ARGS,
    log:
        **logs("05_EVM", "12_evm_stats"),
    shell:
        SH(PY("annotation_stats.py") + """ --set EVM={input.evm} --out {output.stats} --out-dist {output.dist}
        """ + PY("source_support.py") + """ --query {input.evm} {params.sources} \\
            --out-support {output.support} --out-combos {output.combos}
        cat {output.support} {output.combos}
        """)


rule final_stats:
    """Final gene set statistics (coding + lncRNA), source agreement, and a structure comparison
    of all annotation sets (BRAKER4, GeMoMa, EviAnn, EVM, refine, final)."""
    input:
        final=FINAL_INTERNAL,
        lnc=EVIANN_LNC,
        sets=[QC_GFF[k] for k, _ in GENEPRED] + [QC_GFF["evm"], QC_GFF["refine"]],
    output:
        stats=O("05_EVM/stats/final_stats.tsv"),
        dist=O("05_EVM/stats/final_distributions.tsv"),
        support=O("05_EVM/stats/final_source_support.tsv"),
        combos=O("05_EVM/stats/final_source_combinations.tsv"),
        compare=O("05_EVM/stats/all_sets_stats.tsv"),
        compare_dist=O("05_EVM/stats/all_sets_distributions.tsv"),
    params:
        script=script("05_EVM", "13_final_stats"),
        sources=SOURCE_ARGS,
        sets=" ".join(f"--set {n}={QC_GFF[k]}" for k, n in
                      GENEPRED + [("evm", "EVM"), ("refine", "refine"), ("final", "final")]),
    log:
        **logs("05_EVM", "13_final_stats"),
    shell:
        SH(PY("annotation_stats.py") + """ --set final={input.final} --set lncRNA={input.lnc} \\
            --out {output.stats} --out-dist {output.dist}
        """ + PY("source_support.py") + """ --query {input.final} {params.sources} \\
            --out-support {output.support} --out-combos {output.combos}
        """ + PY("annotation_stats.py") + """ {params.sets} --out {output.compare} --out-dist {output.compare_dist}
        cat {output.stats} {output.support} {output.compare}
        """)
