"""
End-to-end, ApJ-ready map of the Perseus Molecular Cloud region containing
NGC 1333 and IC 348. The script:
  • Resolves target coordinates from SIMBAD
  • Downloads a ~5° FITS background (DSS2 Red by default) via SkyView
  • Displays with Astropy WCSAxes (proper RA/Dec, ticks, and grid)
  • Annotates NGC 1333 and IC 348
  • Adds a compass (N/E) and a degree→pc scale bar (assuming a distance)
  • (Optional) Overlays contours from a user-provided FITS, reprojection handled

Dependencies (all pip-installable):
  astropy, matplotlib, numpy, astroquery, reproject (optional if using contours)

Notes:
  • Change SURVEY to e.g. "WISE 3.4" (W1), "DSS2 Blue", "2MASS-K", etc.
  • Adjust FOV_DEG (3–5 recommended to mirror Jørgensen+06 Fig. 3)
  • If providing contours, set CONTOUR_FITS to a FITS file path
  • For ApJ figures, save as PDF/EPS at high DPI
"""

from __future__ import annotations

import os
from io import BytesIO
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval
from astropy.wcs import WCS
from astroquery.simbad import Simbad
# External helpers
from astroquery.skyview import SkyView

# Optional for reprojection of contours
try:
    from reproject import reproject_interp
except Exception:
    reproject_interp = None

# -------------------------------
# User-configurable parameters
# -------------------------------
TARGETS = ["NGC 1333", "IC 348"]
SURVEY = "DSS2 Red"                 # Examples: "WISE 3.4", "DSS2 Red", "DSS2 Blue", "2MASS-K"
FOV_DEG = 5.0                        # Field of view in degrees (width=height)
PIXELS = 900                        # Output image size (square), adjust for desired resolution
DISTANCE_PC = 300.0                  # Assumed distance to Perseus (used for scale bar)

# Contour overlay (optional)
CONTOUR_FITS: Optional[str] = None   # e.g., "/path/to/CO_map.fits" or None to skip
CONTOUR_LEVELS: Optional[Sequence[float]] = None  # e.g., [5, 10, 20, 40] (in map units)

# Output
OUTFIG = "perseus_map_ngc1333_ic348.pdf"

# -------------------------------
# Utilities
# -------------------------------

def resolve_targets(names: Sequence[str]) -> list[SkyCoord]:
    """Resolve a list of target names via SIMBAD (J2000 coordinates)."""
    custom = Simbad()
    custom.add_votable_fields("ra(d)", "dec(d)")
    coords = []
    for name in names:
        result = custom.query_object(name)
        if result is None:
            raise ValueError(f"Could not resolve target: {name}")
        ra_deg = float(result["RA_d"][0])
        dec_deg = float(result["DEC_d"][0])
        coords.append(SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="fk5"))
    return coords


def skyview_cutout(center: SkyCoord, fov_deg: float, survey: str, pixels: int) -> tuple[np.ndarray, WCS]:
    """Fetch a FITS cutout from SkyView and return (data, WCS)."""
    imgs = SkyView.get_images(position=center.to_string("hmsdms"),
                              survey=[survey],
                              coordinates="J2000",
                              width=fov_deg*u.deg,
                              height=fov_deg*u.deg,
                              pixels=pixels)
    if not imgs:
        raise RuntimeError(f"SkyView returned no images for {survey}")
    hdu = imgs[0][0]
    data = hdu.data.astype(float)
    wcs = WCS(hdu.header)
    return data, wcs


def hips2fits_cutout(center, fov_deg=5.0, pixels=1800,
                     hips="DSS2/red", projection="TAN"):
    url = "https://alasky.u-strasbg.fr/hips-image-services/hips2fits"
    params = {
        "hips": hips,
        "width": pixels,
        "height": pixels,
        "fov": fov_deg,
        "projection": projection,
        "coordsys": "icrs",
        "ra": center.ra.deg,
        "dec": center.dec.deg,
        "format": "fits"
    }
    r = requests.get(url, params=params)
    r.raise_for_status()
    hdul = fits.open(BytesIO(r.content))
    return hdul[0].data, WCS(hdul[0].header)



def midpoint(coords: Sequence[SkyCoord]) -> SkyCoord:
    """Compute an approximate midpoint SkyCoord for a small group of sources."""
    ra = u.Quantity([c.ra.wrap_at(180*u.deg).to(u.deg).value for c in coords])
    dec = u.Quantity([c.dec.to(u.deg).value for c in coords])
    ra_mid = np.mean(ra)
    dec_mid = np.mean(dec)
    return SkyCoord(ra=ra_mid*u.deg, dec=dec_mid*u.deg, frame="fk5")


def draw_compass(ax, pos: SkyCoord, length_deg: float = 0.5):
    """Draw N and E arrows at sky position `pos` with WCS-aware transforms."""
    tr = ax.get_transform("fk5")
    # North arrow: +Dec
    ax.annotate("", xy=(pos.ra.deg, pos.dec.deg + length_deg), xytext=(pos.ra.deg, pos.dec.deg),
                xycoords=tr, textcoords=tr,
                arrowprops=dict(arrowstyle="-|>", lw=1.5))
    ax.text(pos.ra.deg, pos.dec.deg + length_deg + 0.05*length_deg, "N",
            transform=tr, ha="center", va="bottom", fontsize=10)
    # East arrow: +RA (eastward)
    ax.annotate("", xy=(pos.ra.deg + length_deg, pos.dec.deg), xytext=(pos.ra.deg, pos.dec.deg),
                xycoords=tr, textcoords=tr,
                arrowprops=dict(arrowstyle="-|>", lw=1.5))
    ax.text(pos.ra.deg + length_deg + 0.05*length_deg, pos.dec.deg, "E",
            transform=tr, ha="left", va="center", fontsize=10)


def draw_scalebar(ax, pos: SkyCoord, length_deg: float, distance_pc: float):
    """Draw a horizontal scale bar of given angular length with a pc label."""
    tr = ax.get_transform("fk5")
    y = pos.dec.deg
    x0 = pos.ra.deg - 0.5*length_deg
    x1 = pos.ra.deg + 0.5*length_deg
    ax.plot([x0, x1], [y, y], transform=tr, lw=2)
    # End caps
    ax.plot([x0, x0], [y-0.01*length_deg, y+0.01*length_deg], transform=tr, lw=2)
    ax.plot([x1, x1], [y-0.01*length_deg, y+0.01*length_deg], transform=tr, lw=2)
    # Label (convert deg to pc)
    pc_per_deg = distance_pc * np.pi/180.0
    label = f"{length_deg:g}°  ≈  {length_deg*pc_per_deg:.1f} pc @ {distance_pc:.0f} pc"
    ax.text(pos.ra.deg, y - 0.08*length_deg, label, transform=tr,
            ha="center", va="top", fontsize=9)


def maybe_overlay_contours(ax, bg_wcs: WCS, contour_fits: Optional[str], levels: Optional[Sequence[float]]):
    if contour_fits is None:
        return
    with fits.open(contour_fits) as hdul:
        cd = hdul[0].data
        cw = WCS(hdul[0].header)
    if reproject_interp is None:
        print("reproject not available; attempting to overlay without reprojection")
        data = cd
    else:
        data, _ = reproject_interp((cd, cw), bg_wcs, shape_out=ax.images[0].get_array().shape)
    if levels is None:
        finite = np.isfinite(data)
        vmin, vmax = np.percentile(data[finite], [70, 99])
        levels = np.linspace(vmin, vmax, 6)
    ax.contour(data, transform=ax.get_transform(bg_wcs), levels=levels, linewidths=0.8)


# -------------------------------
# Main
# -------------------------------
if __name__ == "__main__":
    # Resolve targets and set up center
    coords = resolve_targets(TARGETS)
    center = midpoint(coords)

    # Fetch background
    # data, w = skyview_cutout(center, FOV_DEG, SURVEY, PIXELS)

    # Example: DSS2/red background
    data, wcs = hips2fits_cutout(center, fov_deg=5.0, pixels=1800)

    # Display with WCSAxes
    fig = plt.figure(figsize=(6.5, 6.5))  # Adjust for journal column width
    ax = fig.add_subplot(111, projection=wcs)

    # Contrast/stretch for nice nebulosity/IR background
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch())

    im = ax.imshow(data, origin="lower", cmap="gray", norm=norm)

    # Coordinates, ticks, grid
    ax.set_xlabel("Right Ascension (J2000)")
    ax.set_ylabel("Declination (J2000)")
    ax.coords.grid(True, color="white", ls=":", lw=0.7, alpha=0.5)
    ax.coords[0].set_axislabel("RA (J2000)")
    ax.coords[1].set_axislabel("Dec (J2000)")
    ax.coords[0].set_major_formatter("hh:mm")
    ax.coords[1].set_major_formatter("dd:mm")

    # Markers and labels for targets
    for name, c in zip(TARGETS, coords):
        ax.scatter(c.ra.deg, c.dec.deg, transform=ax.get_transform("fk5"),
                   s=35, facecolor="none", edgecolor="yellow", linewidth=1.2)
        ax.text(c.ra.deg, c.dec.deg + 0.08, name, transform=ax.get_transform("fk5"),
                color="yellow", fontsize=9, ha="center", va="bottom")

    # Compass and scale bar in the lower-right corner
    # Choose positions a bit inside the plot edges
    dec_margin = 0.12 * FOV_DEG
    ra_margin = 0.12 * FOV_DEG
    compass_pos = SkyCoord(ra=center.ra - (0.5*FOV_DEG - ra_margin)*u.deg,
                           dec=center.dec - (0.5*FOV_DEG - dec_margin)*u.deg)
    draw_compass(ax, compass_pos, length_deg=0.3)

    scalebar_pos = SkyCoord(ra=center.ra + (0.15*FOV_DEG)*u.deg,
                            dec=center.dec - (0.5*FOV_DEG - dec_margin)*u.deg)
    draw_scalebar(ax, scalebar_pos, length_deg=1.0, distance_pc=DISTANCE_PC)

    # Optional contours
    maybe_overlay_contours(ax, wcs, CONTOUR_FITS, CONTOUR_LEVELS)

    # Title (optional)
    ax.set_title(f"Perseus Molecular Cloud — {SURVEY} ({FOV_DEG}°)")

    # Tight layout & save
    plt.tight_layout()
    fig.savefig(OUTFIG, dpi=300, bbox_inches="tight")
    print(f"Saved figure → {OUTFIG}")
