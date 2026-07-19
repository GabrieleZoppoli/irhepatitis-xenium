"""Centralised paths for the irHepatitis publication pipeline (Python).

Override at runtime by setting environment variables:
    IRHEP_BUNDLE_ROOT   — bundle root (default: detected from this file's location)
    IRHEP_RAW_DIR       — raw Xenium data dir
    IRHEP_PROCESSED_DIR — processed data dir
"""
import os
from pathlib import Path

_THIS_FILE = Path(__file__).resolve()
BUNDLE_ROOT = Path(os.environ.get("IRHEP_BUNDLE_ROOT",
                                  _THIS_FILE.parent.parent))
RAW_DIR = Path(os.environ.get("IRHEP_RAW_DIR", BUNDLE_ROOT / "data" / "raw"))
PROCESSED_DIR = Path(os.environ.get("IRHEP_PROCESSED_DIR",
                                    BUNDLE_ROOT / "data" / "processed"))
FIGURES_DIR = BUNDLE_ROOT / "figures"
TABLES_DIR = BUNDLE_ROOT / "tables"
SCRIPTS_DIR = BUNDLE_ROOT / "scripts"

SAMPLES = {
    "irHep": ["18321_a", "24774_d"],
    "AIH":   ["54023_b", "56784_c"],
}

SAMPLE_SHIFTS = {
    "18321_a": (0, 0),
    "24774_d": (12000, 0),
    "54023_b": (0, -8000),
    "56784_c": (12000, -8000),
}


def sample_raw_dir(sample_id: str) -> Path:
    """Return the raw Xenium directory path for a sample id."""
    return RAW_DIR / f"0027420_{sample_id}"


if not BUNDLE_ROOT.exists():
    raise RuntimeError(
        f"BUNDLE_ROOT not found: {BUNDLE_ROOT} — set IRHEP_BUNDLE_ROOT"
    )
if not RAW_DIR.exists():
    import warnings
    warnings.warn(
        f"RAW_DIR not found: {RAW_DIR} — populate or set IRHEP_RAW_DIR"
    )
