"""The English config / profile files must hold exactly the same settings as the Chinese ones."""
import os

import pytest
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAIRS = [("config/config.example.yaml", "config/config.example.en.yaml"),
         ("profiles/local/config.yaml", "profiles/local/config.en.yaml"),
         ("profiles/slurm/config.yaml", "profiles/slurm/config.en.yaml")]


@pytest.mark.parametrize("zh,en", PAIRS)
def test_english_version_has_same_settings(zh, en):
    a = yaml.safe_load(open(os.path.join(ROOT, zh)))
    b = yaml.safe_load(open(os.path.join(ROOT, en)))
    assert a == b, f"{zh} and {en} differ: keep both in sync"


DOCS = [("README.md", "README.en.md"), ("docs/deployment_guide.md", "docs/deployment_guide.en.md"),
        ("docs/example_slurm.md", "docs/example_slurm.en.md")]


def structure(path):
    """(heading levels, number of table rows per table, number of code blocks): language independent."""
    lines = open(os.path.join(ROOT, path)).read().splitlines()
    heads = [l.split(" ")[0] for l in lines if l.startswith("#") and not l.startswith("#!")]
    tables, n = [], 0
    for l in lines + [""]:
        if l.startswith("|"):
            n += 1
        elif n:
            tables.append(n)
            n = 0
    return heads, tables, sum(l.startswith("```") for l in lines)


@pytest.mark.parametrize("zh,en", DOCS)
def test_english_document_has_same_structure(zh, en):
    assert structure(zh) == structure(en), f"{zh} and {en} differ in sections / tables: keep both in sync"
