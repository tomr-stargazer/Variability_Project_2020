"""
Makes maps of "number of observations in each band per tile".

I'm borrowing a bunch of code from wuvars-perseus/preliminary_figures.py.

First goal - purge this of atpy dependencies.

"""

import os

import numpy as np
import matplotlib.pyplot as plt

import astropy.table
import pandas as pd


def filter_by_tile(photometry_data, ds):
    """
    Splits the input data tables into 16 tiles, avoiding overlap regions.

    Uses the spatial extent of the `maxvars` spreadsheet to figure out
    where the bounds of the tiles are.

    Only returns data pertaining to our 1202 variables.

    Returns
    -------
    tile_tables : list of atpy.Table
        A list of the variable star photometry tables 
        corresponding to each tile.
    tile_spreadsheets : list of atpy.Table
        A list of the variable star variability spreadsheets
        corresponding to each table.

    """

    vp = photometry_data

    max_ra = ds["median"]["RA"].max()
    min_ra = ds["median"]["RA"].min()
    max_dec = ds["median"]["DEC"].max()
    min_dec = ds["median"]["DEC"].min()

    tile_size_ra = (max_ra - min_ra) / 4
    tile_size_dec = (max_dec - min_dec) / 4

    tile_tables = []
    tile_spreadsheets = []

    for ra_i in range(4):

        for dec_j in range(4):

            tile_min_ra = min_ra + tile_size_ra * ra_i + tile_size_ra / 12
            tile_max_ra = min_ra + tile_size_ra * (ra_i + 1) - tile_size_ra / 12

            tile_min_dec = min_dec + tile_size_dec * dec_j + tile_size_dec / 12
            tile_max_dec = min_dec + tile_size_dec * (dec_j + 1) - tile_size_dec / 12

            tile_photometry = vp[
                (vp["RA"] > tile_min_ra)
                & (vp["RA"] < tile_max_ra)
                & (vp["DEC"] > tile_min_dec)
                & (vp["DEC"] < tile_max_dec)
            ]

            tile_spreadsheet = ds[
                (ds["median"]["RA"] > tile_min_ra)
                & (ds["median"]["RA"] < tile_max_ra)
                & (ds["median"]["DEC"] > tile_min_dec)
                & (ds["median"]["DEC"] < tile_max_dec)
            ]

            tile_tables.append(tile_photometry)
            tile_spreadsheets.append(tile_spreadsheet)

            pdb.set_trace()

    return tile_tables, tile_spreadsheets


def f_observing_map(dat, ds, min_nights=50):
    """
    Makes a map of observations and how many J, H, K band datapoints
    each tile has.

    """

    maxvars = ds[
        (ds["variability"]["Stetson_JHK"] > 0.5)
        & (
            (ds["count"]["N_J"] >= min_nights)
            | (ds["count"]["N_K"] >= min_nights)
            | (ds["count"]["N_H"] >= min_nights)
        )
    ]

    max_ra = maxvars["median"]["RA"].max()
    min_ra = maxvars["median"]["RA"].min()
    max_dec = maxvars["median"]["DEC"].max()
    min_dec = maxvars["median"]["DEC"].min()

    tile_size_ra = (max_ra - min_ra) / 4
    tile_size_dec = (max_dec - min_dec) / 4

    tile_spreadsheets = filter_by_tile(dat, maxvars)[1]

    fig = plt.figure(figsize=(6, 6))

    ij_list = [(x, y) for x in range(4) for y in range(4)]

    text_params = {"horizontalalignment": "center", "verticalalignment": "center"}

    # print the nice text and stuff
    for k, ij, tile_spreadsheet in zip(
        list(range(len(tile_spreadsheets))), ij_list, tile_spreadsheets
    ):

        ra_i, dec_j = ij

        tile_ra = min_ra + tile_size_ra * ra_i + tile_size_ra / 2
        tile_dec = min_dec + tile_size_dec * dec_j + tile_size_dec / 2

        plt.text(
            np.degrees(tile_ra),
            np.degrees(tile_dec + tile_size_dec / 4),
            "Tile #%d" % (k + 1),
            **text_params,
        )
        plt.text(
            np.degrees(tile_ra),
            np.degrees(tile_dec + tile_size_dec / 12),
            "J: %3d" % tile_spreadsheet["count"]["N_J"].max(),
            color="b",
            **text_params,
        )
        plt.text(
            np.degrees(tile_ra),
            np.degrees(tile_dec - tile_size_dec / 12),
            "H: %3d" % tile_spreadsheet["count"]["N_H"].max(),
            color="g",
            **text_params,
        )
        plt.text(
            np.degrees(tile_ra),
            np.degrees(tile_dec - tile_size_dec / 4),
            "K: %3d" % tile_spreadsheet["count"]["N_K"].max(),
            color="r",
            **text_params,
        )

    physical_tile_size_ra = (max_ra - min_ra) / 3.88
    physical_tile_size_dec = (max_dec - min_dec) / 3.88

    # make the overlapping rectangles
    for k, ij in zip(list(range(len(ij_list))), ij_list):

        ra_i, dec_j = ij

        southeast_corner_ra = min_ra + 0.94 * physical_tile_size_ra * ra_i
        southeast_corner_dec = min_dec + 0.94 * physical_tile_size_dec * dec_j

        if ij == (0, 0) or ij == (0, 2) or ij == (2, 0) or ij == (2, 2):
            rectangle_params = {"color": "0.85", "zorder": -10}
        else:
            rectangle_params = {"fill": False}

        plt.gca().add_patch(
            plt.Rectangle(
                (np.degrees(southeast_corner_ra), np.degrees(southeast_corner_dec)),
                np.degrees(physical_tile_size_ra),
                np.degrees(physical_tile_size_dec),
                ec="k",
                **rectangle_params,
            )
        )

    northeast_corner = (
        np.degrees(maxvars["median"]["RA"].max() + 0.001),
        np.degrees(maxvars["median"]["DEC"].max() + 0.001),
    )

    southwest_corner = (
        np.degrees(maxvars["median"]["RA"].min() - 0.001),
        np.degrees(maxvars["median"]["DEC"].min() - 0.001),
    )

    plt.xlim(northeast_corner[0], southwest_corner[0])
    plt.ylim(southwest_corner[1], northeast_corner[1])

    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")

    return fig


if __name__ == "__main__":
    """ Here we're going to use wserv7 as a test case. """

    wserv = 7

    raw_root = (
        "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/Raw_Downloads"
    )
    w7_data_raw_p = os.path.join(raw_root, f"WSERV{str(wserv)}.fits.gz",)

    clean_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"
    w7_data_clean_p = os.path.join(
        clean_root,
        f"wserv{str(wserv)}",
        f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5.h5",
    )

    spreadsheet_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/analysis_artifacts"
    w7_spread_raw_p = os.path.join(
        spreadsheet_root,
        f"wserv{str(wserv)}",
        f"WSERV{str(wserv)}_uncleaned_summary_spreadsheet.h5",
    )

    w7_spread_clean_p = os.path.join(
        spreadsheet_root,
        f"wserv{str(wserv)}",
        f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_summary_spreadsheet.h5",
    )

    w7_data_raw = astropy.table.Table.read(w7_data_raw_p)
    w7_data_clean = astropy.table.Table.read(w7_data_clean_p)

    w7_spread_raw = pd.read_hdf(w7_spread_raw_p, key="table")
    w7_spread_clean = pd.read_hdf(w7_spread_clean_p, key="table")

    fig1 = f_observing_map(dat=w7_data_raw, ds=w7_spread_raw)
    fig2 = f_observing_map(dat=w7_data_clean, ds=w7_spread_clean)
