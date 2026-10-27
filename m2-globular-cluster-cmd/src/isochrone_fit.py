"""Distance and metallicity of NGC 7089 (M2) from the calibrated V vs (V-R) CMD.

What the data can and cannot constrain
--------------------------------------
Our photometry reaches V ~ 18, about a magnitude above the main-sequence
turn-off (V ~ 19.3), so the AGE cannot be measured from these data; we adopt
log(age) = 10.10 (12.6 Gyr) for the isochrones. The horizontal branch (HB)
and the red-giant branch (RGB) are well populated and give:

1. DISTANCE from the HB level (a standard candle):
   (m-M)_V = V_HB - M_V(HB),  M_V(HB) = 0.16 [Fe/H] + 0.84   (Harris 1996, 2010 ed.)
   then subtract A_V = 3.1 E(B-V), with E(B-V) = 0.06 (Harris 2010).

2. METALLICITY from the colour of the RGB: with the distance fixed, compare
   the RGB to PARSEC isochrones of different [M/H] (converted from Gaia bands
   to V and R with the same Evans et al. 2018 relations as the calibration)
   and find the [M/H] whose RGB matches the observed median colour.
   PARSEC uses scaled-solar mixtures; for an alpha-enhanced globular cluster
   ([alpha/Fe] = +0.3) [Fe/H] = [M/H] - 0.21 (Salaris et al. 1993).

Uncertainties: bootstrap over stars, plus the photometric zero-point errors
from calibrate.py propagated into the HB level and the RGB colour.
"""
import glob
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from calibrate import gaia_to_VR

warnings.filterwarnings("ignore")
from config import ISOCHRONES, OUT
rng = np.random.default_rng(7089)

EBV = 0.06                  # Harris (2010)
AV = 3.1 * EBV
EVR = 0.62 * EBV            # Cousins R, R_V = 3.1 extinction law
LOG_AGE = 10.10             # 12.6 Gyr, adopted (turn-off not reached)
ALPHA_CORR = 0.21           # [M/H] - [Fe/H] for [alpha/Fe] = +0.3 (Salaris et al. 1993)
MH_RANGE = (-2.3, -1.0)     # isochrone metallicities to consider

RGB_BOX = dict(V=(13.2, 17.8), VR=(0.38, 1.0))
HB_BOX = dict(V=(15.6, 16.6), VR=(-0.15, 0.25))   # flat (blue) part of the HB


def load_isochrone(path, log_age=LOG_AGE):
    cols = ["Zini", "MH", "logAge", "Mini", "int_IMF", "Mass", "logL", "logTe", "logg",
            "label", "mbol", "G", "BP", "RP"]
    d = pd.DataFrame(np.loadtxt(path, comments="#", usecols=range(14)), columns=cols)
    d = d[np.isclose(d.logAge, log_age, atol=1e-3)].copy()
    d["V"], d["R"] = gaia_to_VR(d.G.values, (d.BP - d.RP).values)
    d["VR"] = d.V - d.R
    return d


def stars():
    s = pd.read_csv(OUT / "ngc7089_VR_calibrated.csv")
    # outside the unresolved core, inside the cluster, good photometry
    return s[(s.r_arcmin > 0.8) & (s.r_arcmin < 5.0) & (s.VR_err < 0.1)]


def in_box(s, box):
    return s[s.V.between(*box["V"]) & s.VR.between(*box["VR"])]


def hb_distance(V_hb, feh):
    dm_v = V_hb - (0.16 * feh + 0.84)
    dm0 = dm_v - AV
    return dm_v, dm0, 10 ** ((dm0 + 5) / 5) / 1000


def rgb_offset(rgb, iso, dm0):
    """Median colour offset (observed - model) of RGB stars at fixed distance."""
    m = iso[iso.label.isin([2, 3])]
    v = (m.V + dm0 + AV).values
    c = (m.VR + EVR).values
    o = np.argsort(v)
    sel = (rgb.V.values > v.min()) & (rgb.V.values < v.max())
    return np.median(rgb.VR.values[sel] - np.interp(rgb.V.values[sel], v[o], c[o]))


def metallicity(rgb, isos, dm0, colour_shift=0.0):
    mhs = np.array(sorted(isos))
    off = np.array([rgb_offset(rgb, isos[m], dm0) for m in mhs]) + colour_shift
    # offset decreases with [M/H]; find the zero crossing by linear interpolation
    return float(np.interp(0.0, off[::-1], mhs[::-1])), mhs, off


def main():
    s = stars()
    rgb, hb = in_box(s, RGB_BOX), in_box(s, HB_BOX)
    cal = json.load(open(OUT / "calibration.json"))
    zpV, zpR = cal["V"]["zp_err"], cal["R"]["zp_err"]
    isos = {}
    for f in sorted(glob.glob(str(ISOCHRONES / "*.dat"))):
        with open(f) as fh:  # read [M/H] from the first data line before loading the whole file
            first = next(l for l in fh if not l.startswith("#"))
        if MH_RANGE[0] <= float(first.split()[1]) <= MH_RANGE[1]:
            iso = load_isochrone(f)
            isos[round(iso.MH.iloc[0], 2)] = iso

    # --- best estimate (iterate [Fe/H] <-> HB distance, converges in 2-3 steps)
    V_hb = float(np.median(hb.V))
    feh = -1.6
    for _ in range(5):
        dm_v, dm0, d = hb_distance(V_hb, feh)
        mh, mhs, offsets = metallicity(rgb, isos, dm0)
        feh = mh - ALPHA_CORR

    # --- uncertainties: bootstrap stars + random zero-point shifts
    boot = []
    for _ in range(1000):
        hb_b = hb.sample(len(hb), replace=True, random_state=rng.integers(1e9))
        rgb_b = rgb.sample(len(rgb), replace=True, random_state=rng.integers(1e9))
        dV, dR = rng.normal(0, zpV), rng.normal(0, zpR)
        vhb = np.median(hb_b.V) + dV
        f = feh
        for _ in range(3):
            _, d0b, db = hb_distance(vhb, f)
            m, _, _ = metallicity(rgb_b.assign(V=rgb_b.V + dV), isos, d0b, colour_shift=dV - dR)
            f = m - ALPHA_CORR
        boot.append((vhb, d0b, db, m, f))
    b = np.array(boot)
    err = lambda col: float(np.std(b[:, col]))

    result = {
        "V_HB": round(V_hb, 3), "V_HB_err": round(err(0), 3), "n_HB": int(len(hb)),
        "(m-M)_V": round(dm_v, 2), "(m-M)_0": round(dm0, 2), "(m-M)_0_err": round(err(1), 2),
        "distance_kpc": round(d, 1), "distance_err_kpc": round(err(2), 1),
        "[M/H]_PARSEC": round(mh, 2), "[M/H]_err": round(err(3), 2),
        "[Fe/H]": round(feh, 2), "[Fe/H]_err": round(err(4), 2), "n_RGB": int(len(rgb)),
        "adopted": {"E(B-V)": EBV, "log_age": LOG_AGE, "age_Gyr": 12.6,
                    "M_V(HB)": "0.16[Fe/H] + 0.84 (Harris 2010)", "alpha_correction": ALPHA_CORR},
        "note": "Statistical + zero-point errors only. Systematics (HB calibration ~0.1 mag, "
                "Gaia->VR transformations ~0.05 mag, extinction) are larger.",
        "literature": {"[Fe/H]": -1.65, "(m-M)_V": 15.50, "distance_kpc": 11.5,
                       "ref": "Harris (1996, 2010 edition)"},
        "rgb_colour_offset_vs_MH": {f"{m:.2f}": round(float(o), 4) for m, o in zip(mhs, offsets)},
    }
    json.dump(result, open(OUT / "isochrone_fit.json", "w"), indent=2)
    print(json.dumps({k: v for k, v in result.items() if k != "rgb_colour_offset_vs_MH"}, indent=2))

    # --- figure: CMD with best-matching isochrone and two bracketing ones
    best_mh = min(isos, key=lambda m: abs(m - mh))
    fig, ax = plt.subplots(figsize=(6.2, 7.2))
    ax.scatter(s.VR, s.V, s=9, c="#2b6cb0", alpha=0.75, edgecolor="none",
               label=f"NGC 7089 stars, HRPO 20-inch ({len(s)})")
    for m, col, lw, lab in [(min(isos), "0.6", 0.9, None), (max(isos), "0.6", 0.9, None),
                            (best_mh, "#e53e3e", 1.8, f"PARSEC [M/H] = {best_mh:.2f}, 12.6 Gyr")]:
        seg = isos[m][isos[m].label.isin([1, 2, 3])]
        ax.plot(seg.VR + EVR, seg.V + dm0 + AV, c=col, lw=lw, label=lab)
        if col == "0.6":
            ax.text(seg.VR.iloc[-1] + EVR + 0.01, seg.V.iloc[-1] + dm0 + AV, f"{m:.2f}", fontsize=7, color="0.45")
    ax.axhline(V_hb, xmin=0.02, xmax=0.42, c="#dd6b20", ls="--", lw=1.2, label=f"HB level V = {V_hb:.2f}")
    ax.set_xlim(-0.35, 1.15)
    ax.set_ylim(18.6, 12.4)
    ax.set_xlabel("V − R [mag]")
    ax.set_ylabel("V [mag]")
    ax.set_title("NGC 7089 (M2): colour–magnitude diagram")
    ax.text(0.97, 0.03,
            f"distance ≈ {d:.1f} ± {err(2):.1f} kpc (HB)\n[Fe/H] ≈ {feh:.2f} ± {err(4):.2f} (RGB colour)\n"
            f"literature: 11.5 kpc, [Fe/H] = −1.65",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5)
    ax.annotate("red-giant branch", (0.72, 13.7), (0.80, 12.9), fontsize=9, arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.annotate("horizontal branch", (0.02, 16.15), (-0.3, 15.1), fontsize=9, arrowprops=dict(arrowstyle="-", lw=0.6))
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "cmd_isochrone.png", dpi=180)


if __name__ == "__main__":
    main()
