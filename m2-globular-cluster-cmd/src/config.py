"""Paths used by all scripts (relative to the project folder)."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = Path(os.environ.get("M2_DATA", ROOT / "data"))
ISOCHRONES = DATA / "isochrones"
GAIA = DATA / "gaia_dr2_ngc7089_2arcmin.tsv"
OUT = ROOT / "results"
COURSE_2022 = ROOT / "course_2022"
