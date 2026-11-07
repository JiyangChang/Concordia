"""Minimal FASTA / GFF helpers shared by the workflow scripts (standard library only)."""
import gzip


def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_fasta(path):
    name, seq = None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def write_fasta(fh, header, seq, width=60):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def parse_attrs(field):
    out = {}
    for part in field.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:  # GTF style: key "value"
            k, v = part.split(" ", 1)
            v = v.strip().strip('"')
        else:
            k, v = part, ""
        out[k.strip()] = v
    return out


def format_attrs(d):
    return ";".join(f"{k}={v}" for k, v in d.items())


def read_gff(path):
    """Yield 9-column lists (comments/FASTA section skipped)."""
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) == 9:
                yield cols
