"""
Makes maps of "number of observations in each band per tile".

I'm borrowing a bunch of code from wuvars-perseus/preliminary_figures.py.

First goal - purge this of atpy dependencies.

"""

import os

import numpy as np
import matplotlib.pyplot as plt

import astropy.table


dropbox_bo_data = os.path.expanduser("~/Desktop/Bo_Tom/NGC1333/WSERV7/DATA/")
dropbox_bo_spread = dropbox_bo_data + "spreadsheet/"

spread = astropy.table.Table.read(
    dropbox_bo_spread + "full_data_errorcorrected_ce_spreadsheet.fits"
)
spread_cleaned = astropy.table.Table.read(
    dropbox_bo_spread + "fdece_graded_clipped0.95_scrubbed0.1_dusted0.5_spread.fits"
)

fulldata = astropy.table.Table.read(
    dropbox_bo_data + "full_data_errorcorrected_ce.fits"
)
fulldata_cleaned = astropy.table.Table.read(
    dropbox_bo_data + "fdece_graded_clipped0.95_scrubbed0.1_dusted0.5.fits"
)


def filter_by_tile(photometry_data, spreadsheet):
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

    max_ra = spreadsheet["RA"].max()
    min_ra = spreadsheet["RA"].min()
    max_dec = spreadsheet["DEC"].max()
    min_dec = spreadsheet["DEC"].min()

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

            tile_spreadsheet = spreadsheet[
                (spreadsheet["RA"] > tile_min_ra)
                & (spreadsheet["RA"] < tile_max_ra)
                & (spreadsheet["DEC"] > tile_min_dec)
                & (spreadsheet["DEC"] < tile_max_dec)
            ]

            tile_tables.append(tile_photometry)
            tile_spreadsheets.append(tile_spreadsheet)

    return tile_tables, tile_spreadsheets


def f_observing_map(data_table=fulldata_cleaned, spreadsheet=spread_cleaned):
    """
    Makes a map of observations and how many J, H, K band datapoints
    each tile has.

    """

    maxvars = spreadsheet[
        (spreadsheet["Stetson"] > 0.5)
        & (
            (spreadsheet["N_j"] >= 50)
            | (spreadsheet["N_k"] >= 50)
            | (spreadsheet["N_h"] >= 50)
        )
    ]

    max_ra = maxvars["RA"].max()
    min_ra = maxvars["RA"].min()
    max_dec = maxvars["DEC"].max()
    min_dec = maxvars["DEC"].min()

    tile_size_ra = (max_ra - min_ra) / 4
    tile_size_dec = (max_dec - min_dec) / 4

    tile_spreadsheets = filter_by_tile(data_table, maxvars)[1]

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
            **text_params
        )
        plt.text(
            np.degrees(tile_ra),
            np.degrees(tile_dec + tile_size_dec / 12),
            "J: %3d" % tile_spreadsheet["N_j"].max(),
            color="b",
            **text_params
        )
        plt.text(
            np.degrees(tile_ra),
            np.degrees(tile_dec - tile_size_dec / 12),
            "H: %3d" % tile_spreadsheet["N_h"].max(),
            color="g",
            **text_params
        )
        plt.text(
            np.degrees(tile_ra),
            np.degrees(tile_dec - tile_size_dec / 4),
            "K: %3d" % tile_spreadsheet["N_k"].max(),
            color="r",
            **text_params
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
                **rectangle_params
            )
        )

    northeast_corner = (
        np.degrees(maxvars["RA"].max() + 0.001),
        np.degrees(maxvars["DEC"].max() + 0.001),
    )

    southwest_corner = (
        np.degrees(maxvars["RA"].min() - 0.001),
        np.degrees(maxvars["DEC"].min() - 0.001),
    )

    plt.xlim(northeast_corner[0], southwest_corner[0])
    plt.ylim(southwest_corner[1], northeast_corner[1])

    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")

    return fig
