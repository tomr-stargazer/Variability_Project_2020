"""
In this script I am coding up maps of the regions we have studied. 

Copied from `map_figure_Perseus_regions.py`.

"""
import os
import warnings
from io import BytesIO

import astropy.constants as c
import astropy.units as u
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import wuvars.analysis.variability_selection as sv
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization import (AsinhStretch, ImageNormalize,
                                   MinMaxInterval, PercentileInterval,
                                   ZScaleInterval, simple_norm)
from astropy.wcs import WCS
from astroquery.simbad import Simbad
from matplotlib.patches import ConnectionPatch
from wuvars.analysis.bd_matching_ic348 import lowmass_ic_joint_matches
from wuvars.analysis.bd_matching_ngc1333 import lowmass_ngc_joint_matches
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.spectral_type_to_number import get_SpT_from_num
from wuvars.analysis.variability_selection_curved import (curve_Stetson, sv_hk,
                                                          sv_jh, sv_jhk, sv_jk)
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.publication_figures.make_figures import figure_export_path

warnings.filterwarnings("ignore")

ngc_match = match_ngc()
ic_match = match_ic()

names = ["ngc", "ic"]
wserv_dict = {"ngc": 7, "ic": 8}
fullname_dict = {"ngc": "NGC 1333", "ic": "IC 348"}
match_dict = {"ngc": ngc_match, "ic": ic_match}
spread_dict = {"ngc": spreadsheet.load_wserv_v2(7), "ic": spreadsheet.load_wserv_v2(8)}
q_dict = {"ngc": quality_classes.load_q(7), "ic": quality_classes.load_q(8)}

make_figs = True

CACHE_DIR = "wise_fits_cache"
os.makedirs(CACHE_DIR, exist_ok=True)


# lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"
# ngc_spread = spreadsheet.load_wserv_v2(7)
# ngc_q = quality_classes.load_q(7)
# inspect_ngc = pd.read_csv(
#     os.path.join(lc_dir, "inspection_ngc.csv"), skipinitialspace=True
# )

# rejected_sources_ngc = inspect_ngc["SOURCEID"][inspect_ngc["exclude?"] == "yes"]
# # inspect_onc['SOURCEID']

# rejected_sources_ngc_ind = np.in1d(
#     lowmass_ngc_joint_matches["SOURCEID"], rejected_sources_ngc
# )
# === IC 348 ===


# lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"
# inspect_ic = pd.read_csv(
#     os.path.join(lc_dir, "inspection_ic.csv"), skipinitialspace=True
# )

# ic_spread = spreadsheet.load_wserv_v2(8)
# ic_q = quality_classes.load_q(8)

# rejected_sources_ic = inspect_ic["SOURCEID"][inspect_ic["exclude?"] == "yes"]
# # inspect_onc['SOURCEID']

# rejected_sources_ic_ind = np.in1d(
#     lowmass_ic_joint_matches["SOURCEID"], rejected_sources_ic
# )


def maps_v4():
    fig = plt.figure(figsize=(5 * 2, 5), dpi=150)
    ax_ngc = fig.add_subplot(122)
    ax_ic = fig.add_subplot(121)

    axes_dict = {"ngc": ax_ngc, "ic": ax_ic}

    for name in names:

        ax_map = axes_dict[name]
        match = match_dict[name]
        spread = spread_dict[name]
        q = q_dict[name]

        ax_map.plot(
            np.degrees(spread["median"]["RA"]),
            np.degrees(spread["median"]["DEC"]),
            "k,",
            alpha=0.25,
            rasterized=True,
        )

        ax_map.plot(
            np.degrees(match.approved["median_RA"]),
            np.degrees(match.approved["median_DEC"]),
            "k.",
            ms=4,
            label=f"Good match (n={len(match.approved['median_RA'])})",
        )
        ax_map.plot(
            np.degrees(match.rejected["median_RA"]),
            np.degrees(match.rejected["median_DEC"]),
            "rx",
            ms=6,
            mew=1.2,
            label=f"Unusable lightcurve (n={len(match.rejected['median_RA'])}))",
        )
        ax_map.invert_xaxis()
        ax_map.set_xlabel("RA (deg)")
        ax_map.set_ylabel("Dec (deg)")

        ax_map.set_xlim(
            np.degrees(spread["median"]["RA"][q.q2]).max(),
            np.degrees(spread["median"]["RA"][q.q2]).min(),
        )
        ax_map.set_ylim(
            np.degrees(spread["median"]["DEC"]).min(),
            np.degrees(spread["median"]["DEC"]).max(),
        )

        ax_map.legend(fontsize=9, loc="upper left")
        ax_map.set_title("NGC 1333", fontsize=18)

    ic_unmatched_coords = np.array(
        [
            (55.99654167, 32.04758333),
            (56.04529167, 32.20408333),
            (56.11066667, 32.13905556),
            (56.17625, 32.20786111),
        ]
    )

    # for pair in ic_unmatched_coords:
    ax_ic.plot(
        ic_unmatched_coords[:, 0],
        ic_unmatched_coords[:, 1],
        color="b",
        linestyle="none",
        marker=(6, 2, 0),
        ms=8,
        mew=1,
        label="No close match (n=4)",
    )

    ax_ic.legend(fontsize=9, loc="upper left")
    ax_ic.set_title("IC 348", fontsize=18)

    ax_ngc.set_aspect(1 / np.cos(np.radians(31)))
    ax_ic.set_aspect(1 / np.cos(np.radians(31)))

    plt.savefig("new_test_IC348_NGC1333_map.png", bbox_inches="tight")
    plt.savefig("new_test_IC348_NGC1333_map.pdf", bbox_inches="tight")

    return fig


def resolve_targets(names):
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


def midpoint(coords):
    ra = u.Quantity([c.ra.wrap_at(180 * u.deg).to(u.deg).value for c in coords])
    dec = u.Quantity([c.dec.to(u.deg).value for c in coords])
    ra_mid = np.mean(ra)
    dec_mid = np.mean(dec)
    return SkyCoord(ra=ra_mid * u.deg, dec=dec_mid * u.deg, frame="fk5")


def hips2fits_cutout(center, fov_deg, pixels, hips, projection="TAN"):
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


def maps_v4_plus_WISEcolor():
    """
    Doing something crazy: bolting the previous maps onto a new figure...

    """

    fig = plt.figure(figsize=(5 * 2, 5 * 2), dpi=150)
    ax_ngc = fig.add_subplot(224)
    ax_ic = fig.add_subplot(223)

    ra_ngc1333, dec_ngc1333 = 52.259, 31.263
    ra_ic348, dec_ic348 = 56.135, 32.158

    FOV_DEG = 7.0  # Field of view in degrees (width=height)
    PIXELS = 1800  # Output image size (square)
    DISTANCE_PC = 300.0  # Assumed distance to Perseus (used for scale bar)

    TARGETS = ["NGC 1333", "IC 348"]

    coords = resolve_targets(TARGETS)
    center = midpoint(coords)

    # Put the new stuff here with the WISE and all that
    rgb_image, wcs = make_wise_rgb_cached(
        center, FOV_DEG, PIXELS, percentile_interval=99.5, asinh_a=0.1
    )

    ax_per = fig.add_subplot(211, projection=wcs)

    ax_per.imshow(rgb_image, origin="lower")
    # ax_per.coords[1].set_limits(30.5, 33.5)

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
    ax_per.set_xlim(x_min_pixel + 60, x_max_pixel - 20)  # hacky
    ax_per.set_ylim(y_min_pixel, y_max_pixel)
    # ax_per.set_xlim(x_min_pixel, x_max_pixel)
    # ax_per.set_ylim(y_min_pixel, y_max_pixel)

    ax_per.set_xlabel("RA (J2000)")
    ax_per.set_ylabel("Dec (J2000)")
    ax_per.coords.grid(True, ls=":", alpha=0.5)
    ax_per.set_title("Perseus Molecular Cloud — WISE RGB (W3/W2/W1)")

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
            transform=ax_per.get_transform("icrs"),
            edgecolor=b["color"],
            facecolor="none",
            lw=2,
            label=b["label"],
            zorder=100
        )
        ax_per.add_patch(rect)

    # outfile = f"{outdir}/WISE_3color.pdf"
    # plt.savefig(outfile, dpi=300, bbox_inches="tight")

    # Keep the below intact -- it is the record of what we have done

    axes_dict = {"ngc": ax_ngc, "ic": ax_ic}

    for name in names:

        ax_map = axes_dict[name]
        match = match_dict[name]
        spread = spread_dict[name]
        q = q_dict[name]

        ax_map.plot(
            np.degrees(spread["median"]["RA"]),
            np.degrees(spread["median"]["DEC"]),
            "k,",
            alpha=0.25,
            rasterized=True,
        )

        ax_map.plot(
            np.degrees(match.approved["median_RA"]),
            np.degrees(match.approved["median_DEC"]),
            "k.",
            ms=4,
            label=f"Good match (n={len(match.approved['median_RA'])})",
        )
        ax_map.plot(
            np.degrees(match.rejected["median_RA"]),
            np.degrees(match.rejected["median_DEC"]),
            "rx",
            ms=6,
            mew=1.2,
            label=f"Unusable lightcurve (n={len(match.rejected['median_RA'])}))",
        )
        ax_map.invert_xaxis()
        ax_map.set_xlabel("RA (deg)")
        ax_map.set_ylabel("Dec (deg)")

        ax_map.set_xlim(
            np.degrees(spread["median"]["RA"][q.q2]).max(),
            np.degrees(spread["median"]["RA"][q.q2]).min(),
        )
        ax_map.set_ylim(
            np.degrees(spread["median"]["DEC"]).min(),
            np.degrees(spread["median"]["DEC"]).max(),
        )

        ax_map.legend(fontsize=9, loc="upper left")
        ax_map.set_title("NGC 1333", fontsize=18)

    ic_unmatched_coords = np.array(
        [
            (55.99654167, 32.04758333),
            (56.04529167, 32.20408333),
            (56.11066667, 32.13905556),
            (56.17625, 32.20786111),
        ]
    )

    # for pair in ic_unmatched_coords:
    ax_ic.plot(
        ic_unmatched_coords[:, 0],
        ic_unmatched_coords[:, 1],
        color="b",
        linestyle="none",
        marker=(6, 2, 0),
        ms=8,
        mew=1,
        label="No close match (n=4)",
    )

    ax_ic.legend(fontsize=9, loc="upper left")
    ax_ic.set_title("IC 348", fontsize=18)

    ax_ngc.set_aspect(1 / np.cos(np.radians(31)))
    ax_ic.set_aspect(1 / np.cos(np.radians(31)))

    ### WE want some lines.

    # Let's break this down into pieces. TO use

    four_ra_dec_pairs = [
        (boxes[1]['ra'][0], boxes[1]['dec'][1]),
        (boxes[1]['ra'][1], boxes[1]['dec'][1]),
        (boxes[0]['ra'][0], boxes[0]['dec'][1]),
        (boxes[0]['ra'][1], boxes[0]['dec'][1]),
    ]

    four_axes = [ax_ic, ax_ic, ax_ngc, ax_ngc]

    for (ra_dec, ax) in zip(four_ra_dec_pairs, four_axes):

        xyA = ax_per.wcs.world_to_pixel_values(*ra_dec)
        coordsA = "data"

        xyB = ra_dec
        coordsB = "data"

        con = ConnectionPatch(
            xyA=xyA, coordsA=coordsA, axesA=ax_per,
            xyB=xyB, coordsB=coordsB, axesB=ax,
            color='0.5', lw=0.5, linestyle='-'
            )
        fig.add_artist(con)

    # fig.plot([0.25, 0.2], [0.75, 0.25], "-", c="0.5", lw=0.5, transform=fig.transFigure)

    # ax_per.transData.transform((x_data, y_data))

    # TODO: rename, move these

    from wuvars.publication_figures.make_figures import figure_export_path


    plt.savefig("XXX_new_test_IC348_NGC1333_map.png", bbox_inches="tight")
    # plt.savefig("XXX_new_test_IC348_NGC1333_map.pdf", bbox_inches="tight")

    fig.savefig(
        os.path.join(figure_export_path, "Figure_2_Map_Perseus.pdf"), bbox_inches="tight"
    )

    return fig


if __name__ == "__main__":
    if True:
        fig = maps_v4_plus_WISEcolor()
    if False:
        fig = maps_v4()
        # legacy_NGC1333_map_figure()
        # legacy_IC348_map_figure()
        # maps_v2()

    if False:
        from wuvars.analysis.bd_matching_v3 import match_ngc, match_ic

        ngc_match = match_ngc()
        ic_match = match_ic()

        print(len(ngc_match.approved))

        fig = plt.figure(figsize=(5 * 2, 5), dpi=150)
        ax_ngc = fig.add_subplot(122)
        ax_ic = fig.add_subplot(121)

        ax_ngc.plot(
            np.degrees(ngc_spread["median"]["RA"]),
            np.degrees(ngc_spread["median"]["DEC"]),
            "k,",
            alpha=0.25,
            rasterized=True,
        )

        ax_ngc.plot(
            np.degrees(ngc_match.approved["median_RA"]),
            np.degrees(ngc_match.approved["median_DEC"]),
            "k.",
            ms=4,
            label=f"Good match (n={len(ngc_match.approved)})",
        )
        ax_ngc.plot(
            np.degrees(ngc_match.rejected["median_RA"]),
            np.degrees(ngc_match.rejected["median_DEC"]),
            "rx",
            ms=6,
            mew=1.2,
            label=f"Unusable lightcurve (n={len(ngc_match.rejected)})",
        )
        ax_ngc.invert_xaxis()
        ax_ngc.set_xlabel("RA (deg)")
        ax_ngc.set_ylabel("Dec (deg)")

        ax_ngc.set_xlim(
            np.degrees(ngc_spread["median"]["RA"][ngc_q.q2]).max(),
            np.degrees(ngc_spread["median"]["RA"][ngc_q.q2]).min(),
        )
        ax_ngc.set_ylim(
            np.degrees(ngc_spread["median"]["DEC"]).min(),
            np.degrees(ngc_spread["median"]["DEC"]).max(),
        )

        ax_ngc.legend(fontsize=9, loc="upper left")
        ax_ngc.set_title("NGC 1333", fontsize=18)

        ax_ic.plot(
            np.degrees(ic_spread["median"]["RA"]),
            np.degrees(ic_spread["median"]["DEC"]),
            "k,",
            alpha=0.25,
            rasterized=True,
        )

        ax_ic.plot(
            np.degrees(ic_match.approved["median_RA"]),
            np.degrees(ic_match.approved["median_DEC"]),
            "k.",
            ms=4,
            label=f"Good match (n={len(ic_match.approved)})",
        )
        ax_ic.plot(
            np.degrees(ic_match.rejected["median_RA"]),
            np.degrees(ic_match.rejected["median_DEC"]),
            "rx",
            ms=6,
            mew=1.2,
            label=f"Unusable lightcurve (n={len(ic_match.rejected)})",
        )
        ax_ic.invert_xaxis()
        ax_ic.set_xlabel("RA (deg)")
        ax_ic.set_ylabel("Dec (deg)")

        ax_ic.set_xlim(
            np.degrees(ic_spread["median"]["RA"][ic_q.q2]).max(),
            np.degrees(ic_spread["median"]["RA"][ic_q.q2]).min(),
        )
        ax_ic.set_ylim(
            np.degrees(ic_spread["median"]["DEC"]).min(),
            np.degrees(ic_spread["median"]["DEC"]).max(),
        )

        # ic_unmatched_coords = np.array(
        #     [
        #         (55.99654167, 32.04758333),
        #         (56.04529167, 32.20408333),
        #         (56.11066667, 32.13905556),
        #         (56.17625, 32.20786111),
        #     ]
        # )

        # for pair in ic_unmatched_coords:
        ax_ic.plot(
            np.degrees(ic_match.unmatched["RA"]),
            np.degrees(ic_match.unmatched["DEC"]),
            color="b",
            linestyle="none",
            marker=(6, 2, 0),
            ms=8,
            mew=1,
            label=f"No match w/in 0.5'' (n={len(ic_match.unmatched)})",
        )

        ax_ic.legend(fontsize=9, loc="upper left")
        ax_ic.set_title("IC 348", fontsize=18)

        ax_ngc.set_aspect(1 / np.cos(np.radians(31)))
        ax_ic.set_aspect(1 / np.cos(np.radians(31)))

        plt.savefig("test_IC348_NGC1333_map.png", bbox_inches="tight")
        plt.savefig("test_IC348_NGC1333_map.pdf", bbox_inches="tight")
