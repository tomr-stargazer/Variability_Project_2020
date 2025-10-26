"""
We're writing a function to make a figure (probably for the Appendix)

which showcases observational uncertainty as a function of magnitude in each band.

It's a classic figure.

I particularly want to know at what magnitude, in each band, the uncertainty crosses
0.05 and 0.02 mag, as well as the 90% completeness limit (somehow I have done this before,
but I forget how...)

"""

import os
import pdb

import astropy.table
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.publication_figures.make_figures import (figure_export_path,
                                                     figureset_export_path)

# data I need: spreadsheets!

spread_dict = {"ngc": spreadsheet.load_wserv_v2(7), "ic": spreadsheet.load_wserv_v2(8)}
q_dict = {"ngc": quality_classes.load_q(7), "ic": quality_classes.load_q(8)}

# This part pulls up the source information from our position-matching scheme.
ngc_match = match_ngc()
ic_match = match_ic()

region_keys = ["ngc", "ic", "both"]
wserv_dict = {"ngc": 7, "ic": 8}

# fmt: off
vlines_dict = {
    "ngc": {
        "J": [17.7, 16.6, 15.4], 
        "H": [17.2, 16.1, 14.7], 
        "K": [16.6, 15.4, 14.2],
        },
}
# fmt: on
vlines_dict["ic"] = vlines_dict["ngc"]
vlines_dict["both"] = vlines_dict["ngc"]

completness_dict = {
    "ngc": {
        "J": 19.1,
        "H": 18.3,
        "K": 17.9,
    },
}
completness_dict["ic"] = completness_dict["ngc"]
completness_dict["both"] = completness_dict["ngc"]


for region_key in region_keys:
    if region_key != "both":
        continue

    try:
        wserv = wserv_dict[region_key]
        ds = spread_dict[region_key]
    except KeyError:
        wserv = "7+8"
        ds = pd.concat([spread_dict["ngc"], spread_dict["ic"]], axis=0)

    q2_all_indices_nonvariable = (
        (ds["count"]["N_J"] > 40)
        & (ds["count"]["N_J"] < 160)
        & (ds["count"]["N_H"] > 40)
        & (ds["count"]["N_H"] < 160)
        & (ds["count"]["N_K"] > 40)
        & (ds["count"]["N_K"] < 160)
        & (ds["median"]["JAPERMAG3"] > 11)
        & (ds["median"]["HAPERMAG3"] > 11)
        & (ds["median"]["KAPERMAG3"] > 11)
        & (ds["count"]["N_J_good"] == ds["count"]["N_J"])
        & (ds["count"]["N_H_good"] == ds["count"]["N_H"])
        & (ds["count"]["N_K_good"] == ds["count"]["N_K"])
        & (ds["median"]["PSTAR"] > 0.75)
        & (ds["variability"]["J_red_chisq"] < 2)
        & (ds["variability"]["H_red_chisq"] < 2)
        & (ds["variability"]["K_red_chisq"] < 2)
    )

    with mpl.rc_context(
        {
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "mathtext.fontset": "stix",  # use STIX fonts (Times-compatible)
        }
    ):

        fig, axes = plt.subplots(3, sharex=True, figsize=(6, 5), dpi=150)

        colors = ["b", "g", "r"]

        axes[0].plot(
            ds["median"]["JAPERMAG3"][q2_all_indices_nonvariable],
            ds["std"]["JAPERMAG3"][q2_all_indices_nonvariable],
            "b,",
            rasterized=True,
        )
        axes[1].plot(
            ds["median"]["HAPERMAG3"][q2_all_indices_nonvariable],
            ds["std"]["HAPERMAG3"][q2_all_indices_nonvariable],
            "g,",
            rasterized=True,
        )
        axes[2].plot(
            ds["median"]["KAPERMAG3"][q2_all_indices_nonvariable],
            ds["std"]["KAPERMAG3"][q2_all_indices_nonvariable],
            "r,",
            rasterized=True,
        )

        for ax, band in zip(axes, ["J", "H", "K"]):

            ax.patch.set_alpha(0)

            ax2 = ax.twinx()
            ax2.set_zorder(-2)
            ax2.hist(
                ds["median"][f"{band}APERMAG3"][q2_all_indices_nonvariable],
                facecolor="0.9",
                edgecolor="k",
                histtype="stepfilled",
                linestyle="-",
                lw=1,
                bins=8*5,
                range=(12, 20),
                alpha=0.5,
                # zorder=-20,
            )
            ax2.tick_params(axis="y", colors="0.5")
            if ax is axes[0]:
                ax2.set_ylabel("N", rotation="horizontal", color="0.5", labelpad=10)

            y_levels = [0.05, 0.02, 0.01]

            for y_level, line in zip(y_levels, vlines_dict[region_key][band]):
                ax.hlines(y_level, xmin=0, xmax=line, color="k", lw=0.5)
                ax.vlines(line, ymin=0, ymax=y_level, color="k", lw=0.5)

            ax.axvline(
                completness_dict[region_key][band], color="k", linestyle="--", lw=0.5
            )

            # ax.axhline(0.02, color='k', alpha=0.5, lw=0.5)
            # ax.axhline(0.01, color='k', alpha=0.5, lw=0.5)

            # ax.set_ylim(0, 0.195)

            ax.set_ylabel(f"$\\sigma_{{{band}}}$", rotation="horizontal", fontsize=14)

            ax.semilogy()
            ax.set_ylim(4e-3, 4e-1)
            ax.set_xlim(12, 20)
            # ax.set_xlim(14, 18.5)
            ax.minorticks_on()
            # ax.tick_params(axis='x', which='minor', bottom=True)

        # axes[0].set_ylabel(axes[0].get_ylabel() + "\n(mag)")

        axes[0].text(12.3, 0.05 * 1.1, "0.05 mag", fontsize=6)
        axes[0].text(12.3, 0.02 * 1.1, "0.02 mag", fontsize=6)
        axes[0].text(12.3, 0.01 * 1.1, "0.01 mag", fontsize=6)

        axes[2].set_xlabel("Median magnitude")
        # axes[0].set_title(f"WSERV{wserv}")
        print(f"WSERV{wserv}")
        plt.subplots_adjust(hspace=0)

        filename = "Figure_0_sigma_v_mag.pdf"
        filepath = os.path.join(figure_export_path, filename)
        print(f"Saving to {filepath}...")
        plt.savefig(filepath)
        # plt.close()
