"""
End-to-end, ApJ-ready map of the Perseus Molecular Cloud region containing
NGC 1333 and IC 348. The script:
  • Resolves target coordinates from SIMBAD
  • Downloads a ~5° FITS background (multi-survey via HiPS → hips2fits)
  • Displays with Astropy WCSAxes (proper RA/Dec, ticks, and grid)
  • Annotates NGC 1333 and IC 348
  • Adds a compass (N/E) and a degree→pc scale bar (assuming a distance)
  • (Optional) Overlays contours from a user-provided FITS, reprojection handled

Dependencies (all pip-installable):
  astropy, matplotlib, numpy, astroquery, requests, reproject (optional if using contours)

Notes:
  • Add/remove surveys in the SURVEYS list (HiPS identifiers, e.g. "DSS2/red", "AllWISE/W1")
  • Adjust FOV_DEG (3–5 recommended to mirror Jørgensen+06 Fig. 3)
  • If providing contours, set CONTOUR_FITS to a FITS file path
  • For ApJ figures, save as PDF/EPS at high DPI
"""

from __future__ import annotations

import os
from io import BytesIO
from typing import Optional, Sequence

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization import (AsinhStretch, ImageNormalize,
                                   MinMaxInterval, PercentileInterval,
                                   ZScaleInterval, simple_norm)
from astropy.wcs import WCS
from astroquery.simbad import Simbad
from astroquery.skyview import SkyView

# Optional for reprojection of contours
try:
    from reproject import reproject_interp
except Exception:
    reproject_interp = None


CACHE_DIR = "wise_fits_cache"
os.makedirs(CACHE_DIR, exist_ok=True)

# -------------------------------
# User-configurable parameters
# -------------------------------


ra_ngc1333, dec_ngc1333 = 52.259, 31.263
ra_ic348, dec_ic348 = 56.135, 32.158


TARGETS = ["NGC 1333", "IC 348"]
SURVEYS = ["DSS2/red", "2MASS/K", "AllWISE/W1", "AllWISE/W3", "SFD"]  # new dust map
FOV_DEG = 7.0  # Field of view in degrees (width=height)
PIXELS = 1800  # Output image size (square)
DISTANCE_PC = 300.0  # Assumed distance to Perseus (used for scale bar)

# Contour overlay (optional)
CONTOUR_FITS: Optional[str] = None
CONTOUR_LEVELS: Optional[Sequence[float]] = None

# Output
outdir = "perseus_maps"
OUTDIR = outdir
os.makedirs(outdir, exist_ok=True)

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
        coords.append(SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="fk5"))
    return coords


def hips2fits_cutout(
    center: SkyCoord, fov_deg: float, pixels: int, hips: str, projection: str = "TAN"
) -> tuple[np.ndarray, WCS]:
    """Fetch a FITS cutout from hips2fits and return (data, WCS)."""
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
        "format": "fits",
    }
    r = requests.get(url, params=params)
    r.raise_for_status()
    hdul = fits.open(BytesIO(r.content))
    return hdul[0].data, WCS(hdul[0].header)


def midpoint(coords: Sequence[SkyCoord]) -> SkyCoord:
    ra = u.Quantity([c.ra.wrap_at(180 * u.deg).to(u.deg).value for c in coords])
    dec = u.Quantity([c.dec.to(u.deg).value for c in coords])
    ra_mid = np.mean(ra)
    dec_mid = np.mean(dec)
    return SkyCoord(ra=ra_mid * u.deg, dec=dec_mid * u.deg, frame="fk5")


def draw_compass(ax, pos: SkyCoord, length_deg: float = 0.5):
    tr = ax.get_transform("fk5")
    ax.annotate(
        "",
        xy=(pos.ra.deg, pos.dec.deg + length_deg),
        xytext=(pos.ra.deg, pos.dec.deg),
        xycoords=tr,
        textcoords=tr,
        arrowprops=dict(arrowstyle="-|>", lw=1.5),
    )
    ax.text(
        pos.ra.deg,
        pos.dec.deg + length_deg + 0.05 * length_deg,
        "N",
        transform=tr,
        ha="center",
        va="bottom",
        fontsize=10,
    )
    ax.annotate(
        "",
        xy=(pos.ra.deg + length_deg, pos.dec.deg),
        xytext=(pos.ra.deg, pos.dec.deg),
        xycoords=tr,
        textcoords=tr,
        arrowprops=dict(arrowstyle="-|>", lw=1.5),
    )
    ax.text(
        pos.ra.deg + length_deg + 0.05 * length_deg,
        pos.dec.deg,
        "E",
        transform=tr,
        ha="left",
        va="center",
        fontsize=10,
    )


def draw_scalebar(ax, pos: SkyCoord, length_deg: float, distance_pc: float):
    tr = ax.get_transform("fk5")
    y = pos.dec.deg
    x0 = pos.ra.deg - 0.5 * length_deg
    x1 = pos.ra.deg + 0.5 * length_deg
    ax.plot([x0, x1], [y, y], transform=tr, lw=2)
    ax.plot(
        [x0, x0], [y - 0.01 * length_deg, y + 0.01 * length_deg], transform=tr, lw=2
    )
    ax.plot(
        [x1, x1], [y - 0.01 * length_deg, y + 0.01 * length_deg], transform=tr, lw=2
    )
    pc_per_deg = distance_pc * np.pi / 180.0
    label = f"{length_deg:g}°  ≈  {length_deg*pc_per_deg:.1f} pc @ {distance_pc:.0f} pc"
    ax.text(
        pos.ra.deg,
        y - 0.08 * length_deg,
        label,
        transform=tr,
        ha="center",
        va="top",
        fontsize=9,
    )


def maybe_overlay_contours(
    ax, bg_wcs: WCS, contour_fits: Optional[str], levels: Optional[Sequence[float]]
):
    if contour_fits is None:
        return
    with fits.open(contour_fits) as hdul:
        cd = hdul[0].data
        cw = WCS(hdul[0].header)
    if reproject_interp is None:
        data = cd
    else:
        data, _ = reproject_interp(
            (cd, cw), bg_wcs, shape_out=ax.images[0].get_array().shape
        )
    if levels is None:
        finite = np.isfinite(data)
        vmin, vmax = np.percentile(data[finite], [70, 99])
        levels = np.linspace(vmin, vmax, 6)
    ax.contour(data, transform=ax.get_transform(bg_wcs), levels=levels, linewidths=0.8)


def sfd_cutout(center: SkyCoord, fov_deg: float, pixels: int):
    """
    Fetch a FITS cutout of the SFD dust map using SkyView.
    Returns (data, WCS).
    """
    result = SkyView.get_images(
        position=center,
        survey="SFD",
        coordinates="J2000",
        pixels=pixels,
        width=fov_deg * u.deg,
        height=fov_deg * u.deg,
    )[0]
    data = result[0].data
    wcs = WCS(result[0].header)
    return data, wcs


def make_wise_rgb(center, fov_deg=5.0, pixels=1800):
    surveys = ["AllWISE/W1", "AllWISE/W2", "AllWISE/W3"]
    bands = []

    for s in surveys:
        print(f"Fetching {s}...")
        data, wcs = hips2fits_cutout(center, fov_deg, pixels, s)

        # percentile interval to clip extreme values
        interval = PercentileInterval(99.5)  # keeps the 0.5%-99.5% range
        vmin, vmax = interval.get_limits(data)
        # normalize to 0..1 and apply asinh stretch
        norm_data = AsinhStretch(a=0.1)((data - vmin) / (vmax - vmin))
        bands.append(norm_data)

    # stack into RGB (R=W3, G=W2, B=W1)
    rgb = np.dstack([bands[2], bands[1], bands[0]])
    rgb = np.clip(rgb, 0, 1)
    return rgb, wcs


def get_wise_band(center, fov_deg, pixels, survey):
    """
    Returns (data, WCS) for a single WISE band.
    Downloads the FITS once, then reuses the local file.
    """
    # Make a safe filename
    safe_name = survey.replace("/", "_") + ".fits"
    filepath = os.path.join(CACHE_DIR, safe_name)

    if os.path.exists(filepath):
        # Load from cache
        with fits.open(filepath) as hdul:
            data = hdul[0].data
            wcs = WCS(hdul[0].header)
        print(f"Loaded {survey} from cache.")
    else:
        # Fetch from HiPS service
        data, wcs = hips2fits_cutout(center, fov_deg, pixels, survey)
        # Save to cache
        hdu = fits.PrimaryHDU(data=data, header=wcs.to_header())
        hdu.writeto(filepath, overwrite=True)
        print(f"Downloaded and cached {survey}.")

    return data, wcs


def make_wise_rgb_cached(
    center, fov_deg=5.0, pixels=1800, percentile_interval=99.5, asinh_a=0.1
):
    surveys = ["AllWISE/W1", "AllWISE/W2", "AllWISE/W3"]
    bands = []

    for s in surveys:
        data, wcs = get_wise_band(center, fov_deg, pixels, s)

        # Percentile + asinh stretch
        interval = PercentileInterval(percentile_interval)
        vmin, vmax = interval.get_limits(data)
        norm_data = AsinhStretch(a=asinh_a)((data - vmin) / (vmax - vmin))
        bands.append(norm_data)

    rgb = np.dstack([bands[2], bands[1], bands[0]])  # R=W3, G=W2, B=W1
    rgb = np.clip(rgb, 0, 1)
    return rgb, wcs


# -------------------------------
# Main
# -------------------------------
if __name__ == "__main__":
    os.makedirs(OUTDIR, exist_ok=True)

    coords = resolve_targets(TARGETS)
    center = midpoint(coords)

    SURVEYS = []
    for hips in SURVEYS:
        try:
            print(f"Fetching survey {hips} ...")
            if hips == "SFD":
                data, w = sfd_cutout(center, FOV_DEG, PIXELS)
            else:
                data, w = hips2fits_cutout(center, FOV_DEG, PIXELS, hips)
        except Exception as e:
            print(f"⚠️ Skipping {hips}: {e}")
            continue

        # make figure
        for mode in ["normal", "inverted"]:
            fig = plt.figure(figsize=(6.5, 6.5))
            ax = plt.subplot(projection=w)

            if mode == "normal":
                cmap = "gray"
            else:
                cmap = "gray_r"  # inverted colormap

            # image
            norm = simple_norm(data, "sqrt", percent=99.5)
            ax.imshow(data, origin="lower", cmap=cmap, norm=norm)

            # coords
            ax.set_xlabel("RA (J2000)")
            ax.set_ylabel("Dec (J2000)")
            ax.coords.grid(True, ls=":", alpha=0.5)

            # annotate targets
            ax.scatter(
                ra_ngc1333,
                dec_ngc1333,
                transform=ax.get_transform("icrs"),
                s=80,
                edgecolor="cyan",
                facecolor="none",
                lw=1.5,
                label="NGC 1333",
            )
            ax.scatter(
                ra_ic348,
                dec_ic348,
                transform=ax.get_transform("icrs"),
                s=80,
                edgecolor="magenta",
                facecolor="none",
                lw=1.5,
                label="IC 348",
            )
            ax.legend(loc="lower right")

            # title
            ax.set_title(f"{hips} ({mode})")

            # save
            safe_name = hips.replace("/", "_")
            outfile = f"{outdir}/{safe_name}_{mode}.pdf"
            plt.savefig(outfile, dpi=300, bbox_inches="tight")
            plt.close(fig)

            print(f"✅ Saved {outfile}")

    # WISE3 grayscale

    if False:

        survey = "AllWISE/W3"
        hips = survey
        mode = "inverted"
        fov_deg = 7.0
        pixels = 1800
        data, wcs = get_wise_band(center, fov_deg, pixels, survey)

        fig = plt.figure(figsize=(6.5, 6.5))
        ax = plt.subplot(projection=wcs)

        cmap = "gray_r"  # inverted colormap

        # image
        norm = simple_norm(data, "sqrt", percent=99.5)
        ax.imshow(data, origin="lower", cmap=cmap, norm=norm)

        # coords
        ax.set_xlabel("RA (J2000)")
        ax.set_ylabel("Dec (J2000)")
        ax.coords.grid(True, ls=":", alpha=0.5)

        # title
        ax.set_title(f"{hips} ({mode})")

        corners = wcs.calc_footprint()

        # Doing a thing with y limits.
        dec_min_world = 30.125
        dec_max_world = 33.25
        ra_min_world, ra_max_world = np.max(corners[0:2, 0]), np.min(corners[2:4, 0])

        # Convert world coordinates to pixel coordinates
        # The origin argument should be 0 for standard Python array indexing
        pix_bl = wcs.world_to_pixel_values(ra_min_world, dec_min_world)
        pix_tr = wcs.world_to_pixel_values(ra_max_world, dec_max_world)

        # Extract pixel coordinates for setting limits
        x_min_pixel, y_min_pixel = pix_bl
        x_max_pixel, y_max_pixel = pix_tr

        # Set the limits using the pixel coordinates
        ax.set_xlim(x_min_pixel + 60, x_max_pixel - 20)  # hacky
        ax.set_ylim(y_min_pixel, y_max_pixel)

        boxes = [
            {
                "ra": (52.72826210855102, 51.685525993509465),
                "dec": (30.84072309236323, 31.73142495998769),
                "color": "cyan",
                "label": "NGC 1333",
            },
            {
                "ra": (56.67586712644869, 55.624367919373086),
                "dec": (31.762148439083123, 32.65672727683106),
                "color": "cyan",
                "label": "IC 348",
            },
        ]
        # Overlay boxes
        for b in boxes:
            # Compute lower-left corner in RA/Dec
            ra0, ra1 = sorted(b["ra"])
            dec0, dec1 = sorted(b["dec"])
            width = ra1 - ra0
            height = dec1 - dec0

            rect = patches.Rectangle(
                (ra0, dec0),
                width,
                height,
                transform=ax.get_transform("icrs"),
                edgecolor=b["color"],
                facecolor="none",
                lw=2,
                label=b["label"],
            )
            ax.add_patch(rect)

        # save
        safe_name = hips.replace("/", "_")
        outfile = f"{outdir}/{safe_name}_{mode}_custom.pdf"
        plt.savefig(outfile, dpi=300, bbox_inches="tight")

        print(f"✅ Saved {outfile}")

    # Three color

    # --- Example plotting ---
    # rgb_image, wcs = make_wise_rgb(center, FOV_DEG, PIXELS)
    rgb_image, wcs = make_wise_rgb_cached(
        center, FOV_DEG, PIXELS, percentile_interval=99.5, asinh_a=0.1
    )

    fig = plt.figure(figsize=(6.5, 6.5))
    ax = plt.subplot(projection=wcs)
    ax.imshow(rgb_image, origin="lower")
    # ax.coords[1].set_limits(30.5, 33.5)

    corners = wcs.calc_footprint()

    # Doing a thing with y limits.
    dec_min_world = 30.125
    dec_max_world = 33.25
    ra_min_world, ra_max_world = np.max(corners[0:2, 0]), np.min(corners[2:4, 0])

    # Convert world coordinates to pixel coordinates
    # The origin argument should be 0 for standard Python array indexing
    pix_bl = wcs.world_to_pixel_values(ra_min_world, dec_min_world)
    pix_tr = wcs.world_to_pixel_values(ra_max_world, dec_max_world)

    # Extract pixel coordinates for setting limits
    x_min_pixel, y_min_pixel = pix_bl
    x_max_pixel, y_max_pixel = pix_tr

    # Set the limits using the pixel coordinates
    ax.set_xlim(x_min_pixel + 60, x_max_pixel - 20)  # hacky
    ax.set_ylim(y_min_pixel, y_max_pixel)
    # ax.set_xlim(x_min_pixel, x_max_pixel)
    # ax.set_ylim(y_min_pixel, y_max_pixel)

    ax.set_xlabel("RA (J2000)")
    ax.set_ylabel("Dec (J2000)")
    ax.coords.grid(True, ls=":", alpha=0.5)
    ax.set_title("Perseus Molecular Cloud — WISE RGB (W3/W2/W1)")

    # Define boxes in decimal degrees
    boxes = [
        {
            "ra": (52.72826210855102, 51.685525993509465),
            "dec": (30.84072309236323, 31.73142495998769),
            "color": "cyan",
            "label": "NGC 1333",
        },
        {
            "ra": (56.67586712644869, 55.624367919373086),
            "dec": (31.762148439083123, 32.65672727683106),
            "color": "cyan",
            "label": "IC 348",
        },
    ]
    # Overlay boxes
    for b in boxes:
        # Compute lower-left corner in RA/Dec
        ra0, ra1 = sorted(b["ra"])
        dec0, dec1 = sorted(b["dec"])
        width = ra1 - ra0
        height = dec1 - dec0

        rect = patches.Rectangle(
            (ra0, dec0),
            width,
            height,
            transform=ax.get_transform("icrs"),
            edgecolor=b["color"],
            facecolor="none",
            lw=2,
            label=b["label"],
        )
        ax.add_patch(rect)

    outfile = f"{outdir}/WISE_3color.pdf"
    plt.savefig(outfile, dpi=300, bbox_inches="tight")

    plt.show()
