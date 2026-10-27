"""Photometric calibration of the V and R PSF photometry against Gaia DR3.

The 2022 version calibrated with a single Gaia star (and, through a typo,
used the Gaia R magnitude to calibrate V). Here every usable star in the
field that is also in Gaia is used:

  1. match V and R detections to each other (1.5 arcsec),
  2. match those to Gaia DR2 (data/gaia_dr2_ngc7089_2arcmin.tsv: VizieR I/345/gaia2, 2 arcmin cone
     around the cluster, downloaded for the course in Nov 2022),
  3. convert Gaia G, BP-RP to Johnson V and Cousins R with the Evans et al.
     (2018, A&A 616, A4) DR2 relations (the same ones used in 2022),
  4. fit a zero point per band, V_std = v_inst + ZP_V (same for R), with
     iterative 3-sigma clipping. A colour term is also fitted as a check.

Output: results/ngc7089_VR_calibrated.csv, results/calibration.json, results/calibration.png
"""
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord

warnings.filterwarnings("ignore")
from config import GAIA, OUT

# Cluster centre (Harris 1996, 2010 edition)
CENTER = SkyCoord(323.36258 * u.deg, -0.82325 * u.deg)


def gaia_to_VR(G, bp_rp):
    """Evans et al. (2018) Table A2 (Gaia DR2), valid for -0.5 < BP-RP < 2.75.
    Scatter of the relations: 0.046 mag (V) and 0.048 mag (R)."""
    x = np.asarray(bp_rp)
    g_minus_v = -0.01760 - 0.006860 * x - 0.1732 * x**2
    g_minus_r = -0.003226 + 0.3833 * x - 0.1346 * x**2
    return G - g_minus_v, G - g_minus_r


def read_gaia(path=GAIA):
    lines = [l for l in open(path) if not l.startswith("#") and l.strip()]
    header = lines[0].rstrip("\n").split(";")
    rows = [l.rstrip("\n").split(";") for l in lines[3:]]  # skip units and dashes lines
    df = pd.DataFrame(rows, columns=header)
    for c in ["RA_ICRS", "DE_ICRS", "Gmag", "e_Gmag", "BPmag", "e_BPmag", "RPmag", "e_RPmag", "BP-RP"]:
        df[c] = pd.to_numeric(df[c].str.strip(), errors="coerce")
    return df.dropna(subset=["RA_ICRS", "DE_ICRS", "Gmag"])


def match(a, b, max_sep=1.5):
    ca = SkyCoord(a.ra.values * u.deg, a.dec.values * u.deg)
    cb = SkyCoord(b.ra.values * u.deg, b.dec.values * u.deg)
    idx, sep, _ = ca.match_to_catalog_sky(cb)
    ok = sep.arcsec < max_sep
    return np.where(ok)[0], idx[ok], sep.arcsec[ok]


def isolation(ra, dec, radius=8.0):
    """Distance (arcsec) to the nearest other detected star."""
    c = SkyCoord(ra * u.deg, dec * u.deg)
    _, sep, _ = c.match_to_catalog_sky(c, nthneighbor=2)
    return sep.arcsec


def clipped_zeropoint(y, sigma=3.0, iters=10):
    """Median zero point with iterative sigma clipping (robust scatter from the MAD)."""
    keep = np.isfinite(y)
    for _ in range(iters):
        zp = np.median(y[keep])
        s = 1.4826 * np.median(np.abs(y[keep] - zp))
        new = np.isfinite(y) & (np.abs(y - zp) < sigma * s)
        if (new == keep).all():
            break
        keep = new
    return zp, s, keep, y - zp


def main():
    v = pd.read_csv(OUT / "psfphot_V.csv")
    r = pd.read_csv(OUT / "psfphot_R.csv")

    iv, ir, _ = match(v, r)
    st = pd.DataFrame({
        "ra": v.ra.values[iv], "dec": v.dec.values[iv],
        "v_inst": v.mag_inst.values[iv], "v_inst_err": v.mag_inst_err.values[iv],
        "r_inst": r.mag_inst.values[ir], "r_inst_err": r.mag_inst_err.values[ir],
    })
    st["vr_inst"] = st.v_inst - st.r_inst
    st["r_arcmin"] = SkyCoord(st.ra.values * u.deg, st.dec.values * u.deg).separation(CENTER).arcmin
    st["nn_arcsec"] = isolation(st.ra.values, st.dec.values)

    g = read_gaia()
    i_s, i_g, sep = match(st, g.rename(columns={"RA_ICRS": "ra", "DE_ICRS": "dec"}), max_sep=1.5)
    cal = st.iloc[i_s].copy()
    gg = g.iloc[i_g]
    cal["V_gaia"], cal["R_gaia"] = gaia_to_VR(gg.Gmag.values, gg["BP-RP"].values)
    cal["bp_rp"] = gg["BP-RP"].values
    good = (
        (cal.r_arcmin > 1.0)            # avoid the crowded core
        & (cal.nn_arcsec > 6)           # isolated in our images (seeing FWHM ~4.7")
        & cal.bp_rp.between(-0.3, 2.5)  # inside the transformation's validity range
        & (cal.v_inst_err < 0.05) & (cal.r_inst_err < 0.05)
    )
    cal = cal[good & np.isfinite(cal.V_gaia)]

    zV, sV, kV, resV = clipped_zeropoint((cal.V_gaia - cal.v_inst).values)
    zR, sR, kR, resR = clipped_zeropoint((cal.R_gaia - cal.r_inst).values)
    nV, nR = int(kV.sum()), int(kR.sum())
    # Colour term check: slope of the residuals against instrumental colour
    kV_ct = np.polyfit(cal.vr_inst.values[kV], resV[kV], 1)[0]
    kR_ct = np.polyfit(cal.vr_inst.values[kR], resR[kR], 1)[0]

    st["V"] = st.v_inst + zV
    st["R"] = st.r_inst + zR
    st["V_err"] = np.hypot(st.v_inst_err, sV / np.sqrt(nV))
    st["R_err"] = np.hypot(st.r_inst_err, sR / np.sqrt(nR))
    st["VR"] = st.V - st.R
    st["VR_err"] = np.hypot(st.V_err, st.R_err)
    st.to_csv(OUT / "ngc7089_VR_calibrated.csv", index=False)

    summary = {
        "gaia_catalogue": "Gaia DR2 (VizieR I/345/gaia2), 2 arcmin cone",
        "n_matched_VR": int(len(st)),
        "n_gaia_matches": int(len(i_s)),
        "n_calibrators": {"V": nV, "R": nR},
        "V": {"zp": round(float(zV), 4), "zp_err": round(float(sV / np.sqrt(nV)), 4),
              "scatter": round(float(sV), 4), "residual_colour_slope": round(float(kV_ct), 3)},
        "R": {"zp": round(float(zR), 4), "zp_err": round(float(sR / np.sqrt(nR)), 4),
              "scatter": round(float(sR), 4), "residual_colour_slope": round(float(kR_ct), 3)},
    }
    json.dump(summary, open(OUT / "calibration.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))

    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    for a, res, keep, band, s in [(ax[0], resV, kV, "V", sV), (ax[1], resR, kR, "R", sR)]:
        x = cal.V_gaia if band == "V" else cal.R_gaia
        a.scatter(x[keep], res[keep], s=10, c="tab:green", label=f"used ({keep.sum()})")
        a.scatter(x[~keep], res[~keep], s=10, c="0.6", marker="x", label="clipped")
        a.axhline(0, c="k", lw=0.8)
        a.set_xlabel(f"{band} from Gaia DR2 [mag]")
        a.set_ylabel(f"{band} residual [mag]")
        a.set_title(f"{band}: scatter {s:.3f} mag")
        a.set_ylim(-0.5, 0.5)
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "calibration.png", dpi=150)


if __name__ == "__main__":
    main()
