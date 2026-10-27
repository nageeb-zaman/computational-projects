"""PSF photometry of the stacked HRPO images of NGC 7089 (M2).

Re-analysis of the ASTR 7361 (Fall 2022) final project.
For each band: subtract a smooth background, build an effective PSF (ePSF)
from bright isolated stars *in that band*, run iterative PSF photometry,
and attach sky coordinates from the plate solution.

Output: results/psfphot_<band>.csv with x, y, RA, Dec, flux, flux_err,
instrumental magnitude and its error.
"""
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.nddata import NDData
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from photutils.background import Background2D, MedianBackground, LocalBackground, MMMBackground
from photutils.detection import DAOStarFinder
from photutils.psf import EPSFBuilder, IterativePSFPhotometry, SourceGrouper, extract_stars

warnings.filterwarnings("ignore")

from config import DATA, OUT


def load(band):
    hdu = fits.open(DATA / f"light_NGC7089_{band}.fit")[1]
    return hdu.data.astype(float), WCS(hdu.header), hdu.header


def subtract_background(img):
    bkg = Background2D(img, box_size=64, filter_size=3,
                       sigma_clip=SigmaClip(sigma=3.0), bkg_estimator=MedianBackground())
    return img - bkg.background, bkg.background_rms_median


def psf_stars(img, rms, fwhm, size=31, n_max=40, min_sep=30, edge=30, core=(512, 512), core_r=180):
    """Bright, isolated, unsaturated stars away from the edges and the crowded core."""
    finder = DAOStarFinder(threshold=30 * rms, fwhm=fwhm)
    src = finder(img).to_pandas()
    xy = src[["x_centroid", "y_centroid"]].to_numpy()
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    np.fill_diagonal(d, np.inf)
    ok = (
        (d.min(axis=1) > min_sep)
        & (src.x_centroid > edge) & (src.x_centroid < img.shape[1] - edge)
        & (src.y_centroid > edge) & (src.y_centroid < img.shape[0] - edge)
        & (np.hypot(src.x_centroid - core[0], src.y_centroid - core[1]) > core_r)
        & (src.peak < np.percentile(src.peak, 95))  # drop the very brightest (possible non-linearity)
    )
    sel = src[ok].sort_values("flux", ascending=False).head(n_max)
    return Table({"x": sel.x_centroid.to_numpy(), "y": sel.y_centroid.to_numpy()})


def run(band, fwhm=7.5, threshold=5.0):
    """fwhm in pixels: measured seeing was ~7.6 px (~4.7 arcsec) in both V and R."""
    img, wcs, hdr = load(band)
    img, rms = subtract_background(img)

    stars = extract_stars(NDData(img), psf_stars(img, rms, fwhm), size=31)
    epsf, _ = EPSFBuilder(oversampling=2, maxiters=10, progress_bar=False)(stars)

    finder = DAOStarFinder(threshold=threshold * rms, fwhm=fwhm, sharplo=0.2, sharphi=1.0, roundlo=-0.8, roundhi=0.8)
    phot = IterativePSFPhotometry(
        epsf, fit_shape=(13, 13), finder=finder, aperture_radius=7,
        grouper=SourceGrouper(min_separation=2.0 * fwhm),
        localbkg_estimator=LocalBackground(15, 25, MMMBackground()),
        maxiters=2,
    )
    # Error map: background noise only (the stack's effective gain is not recorded);
    # the calibration scatter later dominates the bright-star error budget anyway.
    res = phot(img, error=np.full(img.shape, rms))
    df = res.to_pandas()
    df = df[(df.flux_fit > 0) & np.isfinite(df.flux_fit)]
    ra, dec = wcs.pixel_to_world_values(df.x_fit, df.y_fit)
    out = pd.DataFrame({
        "x": df.x_fit, "y": df.y_fit, "ra": ra, "dec": dec,
        "flux": df.flux_fit, "flux_err": df.flux_err,
        "qfit": df.get("qfit", np.nan),
    })
    out["mag_inst"] = -2.5 * np.log10(out.flux)
    out["mag_inst_err"] = 1.0857 * out.flux_err / out.flux
    OUT.mkdir(exist_ok=True)
    out.to_csv(OUT / f"psfphot_{band}.csv", index=False)
    print(f"{band}: background rms {rms:.1f} ADU, {len(stars)} PSF stars, {len(out)} stars measured")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("bands", nargs="*", default=["V", "R"])
    for b in ap.parse_args().bands:
        run(b)
