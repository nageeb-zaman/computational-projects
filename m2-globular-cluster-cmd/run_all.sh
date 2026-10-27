#!/usr/bin/env bash
# Full analysis: PSF photometry -> Gaia calibration -> HB distance + RGB metallicity -> figures
set -e
cd "$(dirname "$0")/src"
python3 photometry.py V R
python3 calibrate.py
python3 isochrone_fit.py
python3 make_comparison.py
