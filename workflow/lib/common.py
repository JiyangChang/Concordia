"""Config loading, validation, sample tables and data-driven decisions (pure Python, no Snakemake).

Everything that decides *which* evidence goes *where* lives here so it can be unit-tested
(tests/test_common.py) independently of Snakemake.
"""

import csv
import os
import re

LONG_PLATFORMS = ("pacbio_isoseq", "pacbio_hifi", "ont_cdna", "ont_drna")
PACBIO = ("pacbio_isoseq", "pacbio_hifi")
ONT = ("ont_cdna", "ont_drna")
STRANDNESS = ("unstranded", "FR", "RF")

MINIMAP2_PRESET = {
    "pacbio_isoseq": "-ax splice:hq -uf",
    "pacbio_hifi": "-ax splice:hq -uf",
    "ont_cdna": "-ax splice",
    "ont_drna": "-ax splice -uf -k14",
}

# GeMoMa ERE strandness keywords
GEMOMA_STRAND = {"unstranded": "FR_UNSTRANDED", "FR": "FR_SECOND_STRAND", "RF": "FR_FIRST_STRAND"}
HISAT2_STRAND = {"unstranded": "", "FR": "--rna-strandness FR", "RF": "--rna-strandness RF"}

# Annotation sources that get BUSCO + OMArk + AGAT QC, and where their qc/ directory lives.
QC_DIRS = {
    "eviann": "02_RNA/EviAnn/qc",
    "braker4": "03_Ab_initio/qc",
    "gemoma": "04_Homology/qc/GeMoMa_merged",   # per-species GeMoMa sets are added in the Snakefile
    "evm": "05_EVM/qc/evm",
    "refine": "05_EVM/qc/refine",
    "final": "05_EVM/qc/final",
}


class ConfigError(ValueError):
    pass


def _abs(path, base):
    """Configured path -> absolute: '~' expanded, relative paths taken from `base` (the working directory)."""
    if path is None or str(path).strip() in ("", "null", "None"):
        return None
    path = os.path.expanduser(str(path))
    return path if os.path.isabs(path) else os.path.normpath(os.path.join(base, path))


def read_table(path):
    """Read a TSV with a header line; skips blank lines and lines starting with '#'."""
    if not path or not os.path.exists(path):
        return []
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.lstrip().startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    return [{k.strip(): (v or "").strip() for k, v in r.items() if k is not None} for r in rows]


# ---------- resources with units ----------
# same conversion as Snakemake's mem resource (humanfriendly sizes, bytes / 1024 / 1024): G = 10^9 bytes = 953.67 MB,
# GiB = 1024 MB; a profile's default-resources mem: "16G" and the config's mem: 16G give the same mem_mb
_MEM_MB = {u: f / 1048576 for u, f in {"m": 1e6, "mb": 1e6, "g": 1e9, "gb": 1e9, "t": 1e12, "tb": 1e12,
                                        "mib": 2**20, "gib": 2**30, "tib": 2**40}.items()}
_MINUTES = {"m": 1, "min": 1, "h": 60, "d": 1440}


def parse_mem_mb(value, default_unit="mb"):
    """'16G', '512M', '1.5T', '64GiB' -> mem_mb as Snakemake computes it (16G -> 15259, 16GiB -> 16384); a bare
    number has `default_unit` ('mb' = MB taken as is, 'gb' = as 'G')."""
    m = re.fullmatch(r"\s*([0-9]*\.?[0-9]+)\s*([A-Za-z]*)\s*", str(value))
    unit = (m.group(2).lower() if m else "") or default_unit
    if not m or unit not in _MEM_MB:
        raise ConfigError(f"memory '{value}': expected a number with unit M, G, T or MiB, GiB, TiB (e.g. 16G)")
    if not m.group(2) and default_unit == "mb":
        return int(float(m.group(1)))                       # mem_mb: 16000 stays 16000
    return max(int(round(float(m.group(1)) * _MEM_MB[unit])), 1)


def parse_minutes(value):
    """'12h', '3d', '1d12h', '90m', or a bare number of minutes -> minutes."""
    v = str(value).strip().lower()
    if re.fullmatch(r"[0-9]+", v):
        return int(v)
    parts = re.findall(r"([0-9]*\.?[0-9]+)\s*(min|d|h|m)", v)
    if not parts or re.sub(r"[0-9.\s]*(min|d|h|m)", "", v):
        raise ConfigError(f"runtime '{value}': expected minutes or a duration such as 12h, 3d, 1d12h")
    return int(round(sum(float(n) * _MINUTES[u] for n, u in parts)))


def normalize_resources(r, name="?"):
    """resources entry -> {threads, mem_mb, runtime}: 'mem' (with unit; bare number = GB) or 'mem_mb' (MB, a unit is
    also accepted), 'runtime' in minutes or with unit."""
    d = dict(r)
    if "mem" in d and "mem_mb" in d:
        raise ConfigError(f"resources.{name}: give either mem or mem_mb, not both")
    if "mem" in d:
        d["mem_mb"] = parse_mem_mb(d.pop("mem"), "gb")
    elif "mem_mb" in d:
        d["mem_mb"] = parse_mem_mb(d["mem_mb"], "mb")
    if "runtime" in d:
        d["runtime"] = parse_minutes(d["runtime"])
    if "threads" in d:
        d["threads"] = int(d["threads"])
    return d


class Project:
    """Validated view of the config plus sample tables."""

    def __init__(self, config, base_dir):
        self.cfg = config
        self.base = base_dir
        self.out = _abs(config.get("outdir", "results"), base_dir)
        self.prefix = config.get("species_prefix") or "GENE"
        self.genome = _abs(config.get("genome"), base_dir)
        self.max_intron = int(config.get("max_intron", 50000))
        rna = config.get("rna", {}) or {}
        self.short = self._load_short(_abs(rna.get("short"), base_dir))
        self.long = self._load_long(_abs(rna.get("long"), base_dir))
        self.rna_allow_none = bool(rna.get("allow_none", False))   # explicit opt-in to annotate without RNA
        self.homologs = self._load_homologs(_abs(config.get("homologs"), base_dir))
        self.orthodb = _abs(config.get("orthodb_proteins"), base_dir)
        sub = config.get("orthodb_subset", {}) or {}
        self.orthodb_levels = [str(x) for x in (sub.get("levels") or [])]
        self.orthodb_levels_tab = _abs(sub.get("levels_tab"), base_dir)
        self.orthodb_level2species_tab = _abs(sub.get("level2species_tab"), base_dir)
        qc = config.get("qc", {}) or {}
        self.busco_lineages = list(qc.get("busco_lineages") or [])
        self.busco_download_path = _abs(qc.get("busco_download_path"), base_dir)
        self.omark_db = _abs(qc.get("omark_db"), base_dir)
        self.reference_gff = _abs(qc.get("reference_gff"), base_dir)
        aq = config.get("assembly_qc", {}) or {}
        self.aq = aq
        mq = aq.get("merqury", {}) or {}
        self.merqury_reads = [_abs(x, base_dir) for x in (mq.get("reads") or [])]
        tel = aq.get("telomere", {}) or {}
        self.telomere_enabled = bool(tel.get("enabled", True))
        # default motif TTTAGGG (plant-type); users may give any other sequence
        self.telomere_motif = (str(tel.get("motif") or "TTTAGGG")).strip().upper()
        fcs = aq.get("fcs_gx", {}) or {}
        self.fcs_db = _abs(fcs.get("db"), base_dir)
        self.fcs_taxid = str(fcs.get("tax_id")) if fcs.get("tax_id") not in (None, "", "null") else None
        self.validate()

    def path(self, value):
        """Any other path setting of the configuration, resolved like the ones above (None when unset)."""
        return _abs(value, self.base)

    # ---------- loading ----------
    def _load_short(self, path):
        rows = []
        for r in read_table(path):
            s = r.get("strandness") or "unstranded"
            if s not in STRANDNESS:
                raise ConfigError(f"rna_short: sample {r.get('sample')} has strandness '{s}', expected one of {STRANDNESS}")
            rows.append({"sample": r["sample"], "r1": _abs(r["r1"], self.base),
                         "r2": _abs(r.get("r2"), self.base), "strandness": s})
        return rows

    def _load_long(self, path):
        rows = []
        for r in read_table(path):
            p = r.get("platform", "")
            if p not in LONG_PLATFORMS:
                raise ConfigError(f"rna_long: sample {r.get('sample')} has platform '{p}', expected one of {LONG_PLATFORMS}")
            rows.append({"sample": r["sample"], "reads": _abs(r["reads"], self.base), "platform": p,
                         # optional: Iso-Seq count table giving each HQ transcript's full-length read count
                         "counts": _abs(r.get("counts") or None, self.base)})
        return rows

    def _load_homologs(self, path):
        rows = []
        for r in read_table(path):
            rows.append({"name": r["name"], "genome": _abs(r["genome"], self.base),
                         "gff3": _abs(r["gff3"], self.base), "protein": _abs(r.get("protein"), self.base)})
        return rows

    # ---------- validation ----------
    def validate(self):
        errs = []
        if not self.genome:
            errs.append("`genome` is required")
        if not self.busco_lineages:
            errs.append("`qc.busco_lineages` must be set by the user (e.g. [eudicotyledons_odb12]); no default is used")
        if not self.short and not self.long and not self.rna_allow_none:
            errs.append("at least one of rna.short / rna.long must list samples "
                        "(set rna.allow_none: true to annotate without RNA-seq: BRAKER4 EP, no EviAnn/UTR/lncRNA)")
        if not self.homologs:
            errs.append("`homologs` must list at least one homologous species")
        if not self.omark_db:
            errs.append("`qc.omark_db` (OMAmer database) is required")
        if self.orthodb_levels and not (self.orthodb and self.orthodb_levels_tab and self.orthodb_level2species_tab):
            errs.append("orthodb_subset.levels needs orthodb_proteins, orthodb_subset.levels_tab and "
                        "orthodb_subset.level2species_tab")
        for h in self.homologs:
            # used in file names, GeMoMa prefixes, miniprot IDs and QC source names
            if not re.match(r"^[A-Za-z0-9_]+$", h["name"]) or h["name"] in ("merged", "ERE", "species"):
                errs.append(f"homolog name '{h['name']}' must contain only [A-Za-z0-9_] and not be "
                            "'merged', 'ERE' or 'species'")
        for tbl, key in ((self.short, "sample"), (self.long, "sample"), (self.homologs, "name")):
            names = [r[key] for r in tbl]
            dup = {n for n in names if names.count(n) > 1}
            if dup:
                errs.append(f"duplicate names: {sorted(dup)}")
        if not set(self.telomere_motif) <= set("ACGT"):
            errs.append("assembly_qc.telomere.motif must contain only A/C/G/T")
        if bool(self.fcs_db) != bool(self.fcs_taxid):
            errs.append("assembly_qc.fcs_gx needs both `db` and `tax_id` (or neither, to skip it)")
        impl = self.evm_cfg.get("implementation", "standard")
        if impl not in ("standard", "rust"):
            errs.append("evm.implementation must be 'standard' or 'rust'")
        for k, r in ((self.cfg.get("resources") or {}).items()):
            try:
                normalize_resources(r or {}, k)
            except ConfigError as e:
                errs.append(f"resources.{k}: {e}")
        if errs:
            raise ConfigError("Config errors:\n  - " + "\n  - ".join(errs))

    # ---------- convenience ----------
    @property
    def evm_cfg(self):
        return self.cfg.get("evm", {}) or {}

    def o(self, *parts):
        return os.path.join(self.out, *parts)

    @property
    def short_names(self):
        return [r["sample"] for r in self.short]

    @property
    def long_names(self):
        return [r["sample"] for r in self.long]

    def short_sample(self, name):
        return next(r for r in self.short if r["sample"] == name)

    def long_sample(self, name):
        return next(r for r in self.long if r["sample"] == name)

    def homolog(self, name):
        return next(r for r in self.homologs if r["name"] == name)

    @property
    def pacbio_names(self):
        return [r["sample"] for r in self.long if r["platform"] in PACBIO]

    @property
    def ont_names(self):
        return [r["sample"] for r in self.long if r["platform"] in ONT]

    @property
    def has_short(self):
        return bool(self.short)

    @property
    def has_long(self):
        return bool(self.long)

    @property
    def has_rna(self):
        return self.has_short or self.has_long

    @property
    def strandness(self):
        """Library strandness used for merged short-read evidence (must be consistent)."""
        vals = {r["strandness"] for r in self.short}
        return vals.pop() if len(vals) == 1 else "unstranded"

    # ---------- assembly QC toggles ----------
    @property
    def run_lai(self):
        return bool(self.aq.get("lai", True))   # computed once EDTA has finished (intact LTRs + RepeatMasker)

    @property
    def run_merqury(self):
        return bool(self.merqury_reads)

    @property
    def run_telomere(self):
        return self.telomere_enabled

    @property
    def run_fcs_gx(self):
        return bool(self.fcs_db and self.fcs_taxid)

    # ---------- decisions ----------
    def merged_bam_groups(self):
        """Which merged BAMs exist: sr, pacbio, ont, lr."""
        g = []
        if self.has_short:
            g.append("sr")
        if self.pacbio_names:
            g.append("pacbio")
        if self.ont_names:
            g.append("ont")
        if self.has_long:
            g.append("lr")
        return g

    def merged_bam_members(self, group):
        if group == "sr":
            return [self.o("02_RNA/align_short", f"{s}.bam") for s in self.short_names]
        if group == "pacbio":
            return [self.o("02_RNA/align_long", f"{s}.bam") for s in self.pacbio_names]
        if group == "ont":
            return [self.o("02_RNA/align_long", f"{s}.bam") for s in self.ont_names]
        if group == "lr":
            return [self.o("02_RNA/align_long", f"{s}.bam") for s in self.long_names]
        raise KeyError(group)

    def braker4_mode(self):
        """BRAKER4 mode from the data: short reads and PacBio long reads are used (dual / ETP / IsoSeq); EP with
        ONT only or without RNA (ONT long reads are never passed to BRAKER4)."""
        sr, pb = self.has_short, bool(self.pacbio_names)
        if sr and pb:
            return "dual"
        if sr:
            return "ETP"
        if pb:
            return "IsoSeq"
        return "EP"

    def braker4_bams(self):
        """(bam_files, isoseq_bam) for BRAKER4 samples.csv; empty string means not used."""
        merged = self.o("02_RNA/merged")
        bam = os.path.join(merged, "sr.bam") if self.has_short else ""
        iso = os.path.join(merged, "pacbio.bam") if self.pacbio_names else ""
        return bam, iso

    def eviann_reads(self):
        """(list of BAM paths, EviAnn tag). Only-ONT is provisionally treated as bam_isoseq."""
        merged = self.o("02_RNA/merged")
        if not self.has_rna:
            return [], "none"
        if self.has_short and self.has_long:
            return [os.path.join(merged, "sr.bam"), os.path.join(merged, "lr.bam")], "bam_mix"
        if self.has_short:
            return [os.path.join(merged, "sr.bam")], "bam"
        return [os.path.join(merged, "lr.bam")], "bam_isoseq"

    def transcript_evidence(self):
        """(stringtie source name, bam group, long-read mode?) for EVM TRANSCRIPT evidence;
        (None, None, False) without RNA."""
        if not self.has_rna:
            return None, None, False
        if self.has_long:
            return "LR_stringtie", "lr", True
        return "SR_stringtie", "sr", False

    def bam_usage(self):
        """Which downstream steps use each merged BAM group (shown in the report)."""
        if not self.has_rna:
            return {}
        use = {g: [] for g in self.merged_bam_groups()}
        bam, iso = self.braker4_bams()
        if bam:
            use["sr"].append("BRAKER4 (bam_files)")
        if iso:
            use["pacbio"].append("BRAKER4 (isoseq_bam)")
        ev_bams, tag = self.eviann_reads()
        for b in ev_bams:
            use[os.path.basename(b)[:-4]].append(f"EviAnn ({tag})")
        tx_src, tx_group, _ = self.transcript_evidence()
        use[tx_group].append(f"StringTie -> EVM TRANSCRIPT ({tx_src})")
        ere = self.gemoma_ere_bam()
        if ere:
            use[os.path.basename(ere)[:-4]].append("GeMoMa ERE")
        if "ont" in use:
            use["ont"].append("merged into lr.bam (not used by BRAKER4)")
        return use

    def gemoma_ere_bam(self):
        merged = self.o("02_RNA/merged")
        if self.has_short:
            return os.path.join(merged, "sr.bam")
        if (self.cfg.get("gemoma", {}) or {}).get("ere_long_reads"):
            return os.path.join(merged, "lr.bam")
        return None

    def evm_overlap(self):
        seg = int(self.evm_cfg.get("segment_size", 500000))
        ov = self.evm_cfg.get("overlap_size")
        if ov:
            ov = int(ov)
        else:
            ov = max(10000, 2 * self.max_intron)
        if ov >= seg:
            ov = seg // 2
        return seg, ov

    def evm_weights(self):
        """Lines of the EVM weights file, keeping only evidence that is actually produced."""
        w = self.evm_cfg.get("weights") or {}
        tx_source = self.transcript_evidence()[0]
        produced = {"BRAKER4", "GeMoMa", "miniprot"} | ({"EviAnn", tx_source} if self.has_rna else set())
        lines = []
        for cls in ("ABINITIO_PREDICTION", "OTHER_PREDICTION", "PROTEIN", "TRANSCRIPT"):
            for src, val in (w.get(cls) or {}).items():
                if src in produced:
                    lines.append(f"{cls}\t{src}\t{val}")
        return lines

    # ---------- EVM weight tuning ----------
    TUNE_DEFAULT_GRID = {"BRAKER4": [1, 2], "GeMoMa": [3, 5], "EviAnn": [5, 10], "miniprot": [2],
                         "transcripts": [5, 10]}

    @property
    def evm_tune(self):
        t = self.evm_cfg.get("tune", {}) or {}
        return {"enabled": bool(t.get("enabled", True)), "region": str(t.get("region", "auto") or "auto"),
                "use_best": bool(t.get("use_best", True)), "grid": t.get("grid") or self.TUNE_DEFAULT_GRID}

    def evm_grid(self):
        """[(combo_id, {source: weight})]: c00 = configured evm.weights, then the grid product.
        'transcripts' in the grid applies to the transcript evidence actually produced
        (LR_stringtie or SR_stringtie). Duplicates of an earlier combination are dropped."""
        import itertools
        tx_source = self.transcript_evidence()[0]
        base = {}
        for line in self.evm_weights():
            _, src, w = line.split("\t")
            base[src] = float(w)
        grid = dict(self.evm_tune["grid"])
        if "transcripts" in grid:
            tx_grid = grid.pop("transcripts")
            if tx_source:
                grid[tx_source] = tx_grid
        keys = [k for k in base if k in grid]
        combos, seen = [("c00", dict(base))], {tuple(sorted(base.items()))}
        for values in itertools.product(*[[float(v) for v in grid[k]] for k in keys]):
            w = dict(base)
            w.update(zip(keys, values))
            key = tuple(sorted(w.items()))
            if key not in seen:
                seen.add(key)
                combos.append((f"c{len(combos):02d}", w))
        return combos

    def evm_weight_lines(self, weights):
        """Weights-file lines for a {source: weight} dict (classes as in evm.weights)."""
        cls = {}
        for line in self.evm_weights():
            c, src, _ = line.split("\t")
            cls[src] = c
        return [f"{cls[s]}\t{s}\t{weights[s]:g}" for s in cls]

    def edta_u_arg(self):
        u = (self.cfg.get("repeat", {}) or {}).get("edta_u")
        return "" if u in (None, "", "null") else f"-u {float(u)}"

    def resources(self, name):
        """threads / mem_mb (MB) / runtime (minutes) of a step: resources.<name> over resources.default."""
        res = self.cfg.get("resources", {}) or {}
        d = normalize_resources(res.get("default") or {"threads": 4, "mem_mb": 8000, "runtime": 600}, "default")
        d.update(normalize_resources(res.get(name) or {}, name))
        return d
