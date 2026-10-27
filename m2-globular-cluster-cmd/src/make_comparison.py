"""Side-by-side figure: the 2022 colour-magnitude diagram and the 2026 re-analysis."""
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from config import COURSE_2022, OUT

OLD = COURSE_2022 / "cmd_V_2022.png"

fig, ax = plt.subplots(1, 2, figsize=(12, 6.2), gridspec_kw={"width_ratios": [1.15, 1]})
ax[0].imshow(mpimg.imread(OLD))
ax[0].set_title("2022 (course project)\nsingle-star calibration, shared PSF", fontsize=11)
ax[1].imshow(mpimg.imread(OUT / "cmd_isochrone.png"))
ax[1].set_title("2026 re-analysis\nper-band ePSF, Gaia calibration on ~70 stars", fontsize=11)
for a in ax:
    a.axis("off")
fig.tight_layout()
fig.savefig(OUT / "before_after.png", dpi=150)
print("wrote", OUT / "before_after.png")
