"""Trim full PARSEC isochrone tables to a single age, so the repository stays small.

The original tables (one file per metallicity, all ages from log t = 5 to 10.45)
were downloaded in 2022 from the CMD web interface (http://stev.oapd.inaf.it/cgi-bin/cmd)
in the Gaia DR2 passbands, about 3 MB each. Only log t = 10.10 (12.6 Gyr) is used here.

usage: python tools/trim_isochrones.py FULL_TABLE.dat [...] --out data/isochrones --logage 10.10
"""
import argparse
from pathlib import Path


def trim(path, out_dir, logage):
    colhdr = next(l for l in open(path) if l.startswith("# Zini"))
    rows = [l for l in open(path)
            if not l.startswith("#") and l.split() and abs(float(l.split()[2]) - logage) < 1e-3]
    mh = float(rows[0].split()[1])
    out = Path(out_dir) / f"parsec_logAge{logage:.2f}_MH{mh:+.2f}.dat"
    with open(out, "w") as fh:
        fh.write("# PARSEC isochrone (CMD web interface, Gaia DR2 passbands), downloaded Nov 2022 for ASTR 7361.\n")
        fh.write(f"# Trimmed to log(age) = {logage:.2f} by tools/trim_isochrones.py. Source file: {Path(path).name}\n")
        fh.write(colhdr)
        fh.writelines(rows)
    return out, len(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--out", default="data/isochrones")
    ap.add_argument("--logage", type=float, default=10.10)
    a = ap.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    for f in a.files:
        print(*trim(f, a.out, a.logage))
