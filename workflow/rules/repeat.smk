# 01_Repeat: EDTA TE library -> RepeatMasker soft-masking.

REP = config.get("repeat", {}) or {}
_cds = REP.get("edta_cds")
if isinstance(_cds, bool):
    raise ValueError("repeat.edta_cds is now a FILE PATH (CDS nucleotide FASTA for EDTA --cds) or null; "
                     f"got {_cds!r}. The homolog CDS are no longer used by default.")
EDTA_CDS = P.path(_cds)


DFAM_DIR = P.path(REP.get("dfam_dir"))


rule repeat_famdb_setup:
    """RepeatMasker >= 4.2 needs a Dfam FamDB (at least the root partition) even with a custom -lib;
    the bioconda package ships none. If famdb.py of the environment finds no data, link the files
    of repeat.dfam_dir into <famdb install dir>/Libraries/famdb (famdb's default location). This
    repeats automatically for every newly created conda environment."""
    output:
        O("01_Repeat/famdb_setup.done"),
    params:
        script=script("01_Repeat", "00_famdb_setup"),
        dfam=DFAM_DIR or "",
    log:
        **logs("01_Repeat", "00_famdb_setup"),
    conda:
        "../envs/repeat.yaml"
    shell:
        SH("""
        FAMDB=$(command -v famdb.py || true)
        if [ -z "$FAMDB" ]; then echo "ERROR: famdb.py not found (RepeatMasker >= 4.2 environment expected)" >&2; exit 1; fi
        FD=$(dirname "$(python3 -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "$FAMDB")")
        # bin/famdb.py may be a wrapper; the installation is the share/famdb-* directory next to it
        [ -f "$FD/famdb_globals.py" ] || FD=$(ls -d "$(dirname "$FAMDB")"/../share/famdb-* 2>/dev/null | head -n 1)
        if python3 "$FD/famdb.py" info > /dev/null 2>&1; then
            echo "FamDB already available in $FD"
        else
            if [ -z "{params.dfam}" ] || ! ls {params.dfam}/*.h5 > /dev/null 2>&1; then
                echo "ERROR: RepeatMasker needs a Dfam FamDB; set repeat.dfam_dir to a directory with dfam*.0.h5 (see docs/deployment_guide.md)" >&2
                exit 1
            fi
            mkdir -p "$FD/Libraries/famdb"
            ln -sf {params.dfam}/*.h5 "$FD/Libraries/famdb/"
            echo "linked $(ls {params.dfam}/*.h5 | tr '\\n' ' ') into $FD/Libraries/famdb"
        fi
        python3 "$FD/famdb.py" info | sed -n '1,8p'      # sed reads everything (head would cause SIGPIPE)
        touch {output}
        """, env="repeat")


rule edta:
    input:
        famdb=rules.repeat_famdb_setup.output,
        genome=GENOME,
        cds=[EDTA_CDS] if EDTA_CDS else [],   # repeat.edta_cds: user-specified CDS FASTA, default none
    output:
        lib=O("01_Repeat/EDTA/genome.fa.mod.EDTA.TElib.fa"),
        anno=O("01_Repeat/EDTA/genome.fa.mod.EDTA.TEanno.gff3"),
        summary=O("01_Repeat/EDTA/genome.fa.mod.EDTA.TEanno.sum"),
    params:
        script=script("01_Repeat", "01_edta"),
        workdir=O("01_Repeat/EDTA"),
        sensitive=REP.get("edta_sensitive", 1),
        cds=lambda w, input: f"--cds {input.cds[0]}" if input.cds else "",
        u=P.edta_u_arg(),  # empty -> EDTA default 1.3e-8 (see config comments)
    threads: res("edta")["threads"]
    resources:
        mem_mb=res("edta")["mem_mb"],
        runtime=res("edta")["runtime"],
    log:
        **logs("01_Repeat", "01_edta"),
    conda:
        "../envs/repeat.yaml"
    shell:
        SH("""
        ln -sf {input.genome} genome.fa
        EDTA.pl --genome genome.fa --species others --step all \\
                --sensitive {params.sensitive} --anno 1 --threads {threads} {params.cds} {params.u}
        """, env="repeat", workdir="{params.workdir}", versions="EDTA.pl -h 2>&1 | grep -i version | head -n 1")


rule repeatmasker:
    input:
        famdb=rules.repeat_famdb_setup.output,
        genome=GENOME,
        lib=rules.edta.output.lib,
    output:
        masked=SOFTMASKED,
        gff=O("01_Repeat/RepeatMasker/genome.fa.out.gff"),
        out=O("01_Repeat/RepeatMasker/genome.fa.out"),
        tbl=O("01_Repeat/RepeatMasker/genome.fa.tbl"),
        divsum=O("01_Repeat/RepeatMasker/genome.divsum"),
    params:
        script=script("01_Repeat", "02_repeatmasker"),
        workdir=O("01_Repeat/RepeatMasker"),
        pa=max(1, res("repeatmasker")["threads"] // 4),  # each RMBlast job uses ~4 cores
        nolow="-nolow" if REP.get("rm_nolow", False) else "",
    threads: res("repeatmasker")["threads"]
    resources:
        mem_mb=res("repeatmasker")["mem_mb"],
        runtime=res("repeatmasker")["runtime"],
    log:
        **logs("01_Repeat", "02_repeatmasker"),
    conda:
        "../envs/repeat.yaml"
    shell:
        SH("""
        ln -sf {input.genome} genome.fa
        RepeatMasker -pa {params.pa} -e rmblast -lib {input.lib} -xsmall {params.nolow} -gff -a -dir . genome.fa
        cp genome.fa.masked {output.masked}
        # bioconda copies util scripts into bin/; they look for RepeatMaskerConfig.pm one level up from
        # themselves, i.e. the env root -> point PERL5LIB at the real RepeatMasker directory
        export PERL5LIB="$(dirname "$(readlink -f "$(command -v RepeatMasker)")")${{PERL5LIB:+:$PERL5LIB}}"
        calcDivergenceFromAlign.pl -s genome.divsum genome.fa.align
        """, env="repeat", workdir="{params.workdir}", versions="RepeatMasker -v")


REPEAT_STATS = ["overview", "te_classes", "te_methods", "library", "landscape", "ltr_age"]


rule repeat_stats:
    """Repeat annotation statistics (EDTA + RepeatMasker) for 07_Report."""
    input:
        te=rules.edta.output.anno,
        lib=rules.edta.output.lib,
        tbl=rules.repeatmasker.output.tbl,
        divsum=rules.repeatmasker.output.divsum,
    output:
        [O(f"01_Repeat/stats/{n}.tsv") for n in REPEAT_STATS],
    params:
        script=script("01_Repeat", "03_repeat_stats"),
        outdir=O("01_Repeat/stats"),
        edta_dir=O("01_Repeat/EDTA"),
        u=(f"{float(REP['edta_u'])} (user)" if REP.get("edta_u") not in (None, "", "null")
           else "1.3e-8 (EDTA default, rice)"),
    log:
        **logs("01_Repeat", "03_repeat_stats"),
    shell:
        SH("""
        # intact-LTR list: EDTA does not document its exact location (same search as the LAI rule)
        PASS=$(ls {params.edta_dir}/*EDTA.raw/LTR/*.pass.list {params.edta_dir}/*EDTA.raw/*.pass.list 2>/dev/null | head -n 1 || true)
        echo "intact LTR list: ${{PASS:-not found}}"
        """ + PY("repeat_stats.py") + """ --te-gff {input.te} --tbl {input.tbl} --lib {input.lib} \\
            --divsum {input.divsum} --pass-list "$PASS" --mutation-rate "{params.u}" --outdir {params.outdir}
        cat {params.outdir}/overview.tsv {params.outdir}/te_classes.tsv
        """)
