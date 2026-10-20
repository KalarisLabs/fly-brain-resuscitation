"""Central path resolution for the fly-brain-resuscitation repository.

Every script resolves its inputs and outputs through three roots, each overridable by an environment variable:

  RAW_ROOT       FlyWire / Codex downloads (connections_princeton*.csv.gz, neurons*.csv.gz, ...)   default: <repo>/data/raw
  RESEARCH_ROOT  working tree for caches, raw run outputs, tables and figures                       default: <repo>/outputs
  FLYSIM_CACHE   FAFB simulator caches (fafb783_*.npz/csv)                                           default: $RESEARCH_ROOT/cache

The layout under RESEARCH_ROOT mirrors the original research tree (sim/, structure/, dynamics/, followup/...),
so pointing RESEARCH_ROOT at an existing research directory reuses its caches and results unchanged.
"""
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("RAW_ROOT", os.path.join(REPO, "data", "raw"))
OUT = os.environ.get("RESEARCH_ROOT", os.path.join(REPO, "outputs"))
CACHE = os.environ.get("FLYSIM_CACHE", os.path.join(OUT, "cache"))


def out(*parts, make=False):
    """Path under RESEARCH_ROOT; make=True creates the parent directory."""
    p = os.path.join(OUT, *parts)
    if make:
        os.makedirs(os.path.dirname(p), exist_ok=True)
    return p


def raw(name):
    """Path of a raw FlyWire download under RAW_ROOT."""
    return os.path.join(RAW, name)
