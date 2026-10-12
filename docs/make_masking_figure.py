"""Draw the README masking walkthrough using a small, deterministic cube.

Run from any directory after installing SpectralCubePlus:
    python /path/to/spectral_cube_plus/docs/make_masking_figure.py
"""
from pathlib import Path
import sys

# Use this checkout even if an older release is installed in the environment.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from matplotlib.patches import Patch, Rectangle
import numpy as np
from scipy import ndimage
import astropy.units as u
from astropy.wcs import WCS
from radio_beam import Beam

from spectral_cube_plus import SpectralCubePlus
from spectral_cube_plus.spectralcubeplus_masks import (
    _beam_geometry, get_expmask, get_prunemask, get_prunemaskvelo,
)


# Constant RMS = 1 K makes the toy intensities equal to signal-to-noise.
# Five channels contain the source; the displayed plane is the middle one.
data = np.full((9, 31, 41), 0.5)
y, x = np.mgrid[:31, :41]
radius = np.hypot(x - 13, y - 14)
data[2:7, radius <= 8] = 3
data[2:7, radius <= 4] = 7
data[2:7, 14, 13] = 0.5  # A valid, below-threshold hole.
data[2:7, (x - 31)**2 + (y - 14)**2 <= 4] = 3  # Detached faint patch.
data[4, 7, 32] = 8  # Bright in only one channel.
data[2:7, 21:23, 31:33] = 7  # Persistent, but too small spatially.

wcs = WCS(naxis=3)
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN", "VRAD"]
wcs.wcs.cunit = ["deg", "deg", "m/s"]
wcs.wcs.cdelt = [-1 / 3600, 1 / 3600, 1000]
wcs.wcs.crpix = [21, 16, 1]
cube = SpectralCubePlus(data * u.K, wcs=wcs,
                        beam=Beam(major=2*u.arcsec, minor=2*u.arcsec))
rms = np.ones(data.shape[1:]) * u.K
options = dict(hthresh=5, lthresh=2, npix=1, radfrac=1, beamarea=2)

# Follow the same stages and helpers as get_expmask, exposing intermediates.
low = data > options["lthresh"]
high = data > options["hthresh"]
seeds = get_prunemaskvelo(high, npix=options["npix"])
grown = ndimage.binary_dilation(seeds, iterations=-1, mask=low)
beam_pixels, kernel = _beam_geometry(cube, options["radfrac"])
closed = ndimage.binary_closing(grown, structure=kernel)
final = get_prunemask(closed, thresh=options["beamarea"] * beam_pixels)
np.testing.assert_array_equal(final, get_expmask(cube, rms, **options))

# Verify that the example really demonstrates every advertised transition.
assert high[4, 7, 32] and not seeds[4, 7, 32]
assert grown[4, 14, 20] and not high[4, 14, 20]
assert not grown[4, 14, 31]  # Faint island has no surviving seed.
assert closed[4, 14, 13] and not grown[4, 14, 13]
assert closed[4, 21, 31] and not final[4, 21, 31]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
orange, blue, purple = "#dc8421", "#247ca8", "#9b57b7"
intensity = LinearSegmentedColormap.from_list(
    "toy_intensity", ["#f5f7f8", "#aeced7", "#16465c"])
fig, axes = plt.subplots(2, 3, figsize=(14, 10))
fig.subplots_adjust(left=.035, right=.985, bottom=.17, top=.85,
                    wspace=.13, hspace=.48)
fig.suptitle("How an expanded emission mask is built", x=.035, y=.975,
             ha="left", fontsize=23, fontweight="bold", color="#203844")
fig.text(.035, .934, "A toy cube: keep reliable bright emission, grow into connected faint emission, then tidy the mask.",
         fontsize=12, color="#52616b")
fig.text(.035, .899, "Same channel in every panel  |  high > 5 RMS  |  low > 2 RMS  |  npix = 1  |  radfrac = 1  |  beamarea = 2",
         fontsize=11, color="#52616b")

titles = ["1  Input data", "2  Find bright seeds", "3  Check neighbouring channels",
          "4  Grow into faint emission", "5  Close small gaps", "6  Remove small spatial patches"]
captions = ["Core = 7 RMS; faint outskirts = 3 RMS.\nHole and background = 0.5 RMS.",
            "Orange samples exceed the high threshold.\nFaint samples are not seeds.",
            "Seeds need high samples on both sides.\nThe one-channel spike loses its seed.",
            "Grow through connected samples > 2 RMS.\nThe detached faint patch stays excluded.",
            "Closing fills the small hole (purple),\neven though it is below 2 RMS.",
            "Cut small areas separately in each channel.\nBlue shows the final inclusion mask."]


def overlay(ax, mask, colour):
    rgba = np.zeros((*mask.shape, 4))
    rgba[mask] = to_rgba(colour)
    ax.imshow(rgba, origin="lower", interpolation="nearest")


for i, ax in enumerate(axes.flat):
    image = ax.imshow(data[4], origin="lower", cmap=intensity, vmin=0, vmax=8,
                      interpolation="nearest", alpha=1 if i == 0 else .35)
    if i in (1, 2):
        overlay(ax, (high if i == 1 else seeds)[4], orange)
    elif i in (3, 4):
        overlay(ax, (grown if i == 3 else closed)[4], blue)
        overlay(ax, seeds[4], orange)
        if i == 4:
            overlay(ax, (closed & ~grown)[4], purple)
    elif i == 5:
        overlay(ax, final[4], blue)
    ax.set_title(titles[i], loc="left", fontsize=13, fontweight="bold", pad=12)
    ax.text(0, -.065, captions[i], transform=ax.transAxes, va="top",
            fontsize=11, linespacing=1.5, color="#52616b")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor("#d7e1e5")

# Labels are confined to the input panel; identical positions aid comparison.
ax = axes[0, 0]
for label, point, text_position in [
    ("bright core", (11, 15), (1, 26)),
    ("faint outskirts", (7, 10), (0, 1)),
    ("hole", (13, 14), (22, 3)),
    ("faint island", (31, 14), (26, 17)),
    ("one-channel spike", (32, 7), (24, 0)),
    ("small bright clump", (32, 22), (21, 27)),
]:
    ax.annotate(label, xy=point, xytext=text_position, fontsize=9,
                color="#203844", arrowprops=dict(arrowstyle="-", color="#617c88"))

# An inset makes the spectral check visible alongside the spatial snapshot.
inset = axes[0, 2].inset_axes([.015, .77, .70, .22])
inset.set_facecolor("white")
inset.set_xlim(-2.9, 3)
inset.set_ylim(-.25, 2.25)
for row, values, label in [(1, [7, 7, 7], "core"), (0, [.5, 8, .5], "spike")]:
    inset.text(-2.7, row+.4, label, va="center", fontsize=9)
    for j, value in enumerate(values):
        inset.add_patch(Rectangle((j-1.1, row), .85, .8,
                                   facecolor=orange if value > 5 else "#e9eff2"))
        inset.text(j-.675, row+.4, f"{value:g}", ha="center", va="center", fontsize=9)
inset.text(-.7, 2, "previous / shown / next channel", fontsize=8)
inset.axis("off")

fig.legend(handles=[Patch(color=orange, label="High-threshold seeds"),
                    Patch(color=blue, label="Included emission"),
                    Patch(color=purple, label="Added by closing")],
           loc="lower left", bbox_to_anchor=(.025, .048), ncol=3,
           frameon=False, fontsize=11)
fig.text(.035, .027, "Seed checking, growth and closing use the 3D cube. Spatial area removal acts on each 2D channel separately.",
         fontsize=11, color="#52616b")
output = Path(__file__).resolve().parents[1] / "docs" / "images"
output.mkdir(parents=True, exist_ok=True)
fig.savefig(output / "expanded_mask_walkthrough.png", dpi=180, facecolor="white")
fig.savefig(output / "expanded_mask_walkthrough.svg", facecolor="white")
plt.close(fig)
print(f"Wrote masking walkthrough to {output}; final mask matches get_expmask.")
