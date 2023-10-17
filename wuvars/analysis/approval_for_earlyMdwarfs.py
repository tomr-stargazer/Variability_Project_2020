"""
Here's the rough cut.

Basically: 
In October 2023 I decided to ~expand~ the number of objects under study in this paper
from just the post-M4.5 objects in NGC 1333 and IC 348 (of which there were ~328 
objects) to also include M0-M4.25 objects (of which there are an additional 221). 

Plus, there are (I think) three objects in IC 348 which had positional match distances
between 0.37 and 0.5 arcsec; I'll try to re-include them here.

The activity being undertaken here:
* I make a bunch of light curves and I save them as graphics files
* I make a spreadsheet that contains information on those objects
* I then, separately, scan through all the lightcurves by eye with the spreadsheet open,
  and note in the "exclude?" column the word "yes" whenever I see a disqualifiable light
  curve.

Update: I have done this. The data live in these two files:
/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3.5/
new_inspection_ngc.csv
new_inspection_ic.csv

I have added the word "yes" to the column "exclude?" wherever I identified a light curve
that was dominated by obvious artifacts (chip edge issues, barely any data, obviously 
unphysical stuff, etc). I tried not to over-use this when I was just looking at a 
not-convincingly-variable Q0 object. I flagged 11 new rejects in IC 348 and one in 
NGC 1333. One of the IC 348 rejects was the 0.49'' match. (Though I did find a 0.50'' 
L0 match which I find convincing. I'll chalk its uncertain position up to its faintness.)

"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.q_string import q_string
from wuvars.analysis.sidsep import sidsep
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.plotting.lightcurve import (ic348_simple_lc_scatter_brokenaxes,
                                        ngc1333_simple_lc_scatter_brokenaxes)

# Let's load up the selection of all the objects we'd like to take a look at.
ngc_match = match_ngc()
ic_match = match_ic()

# Let's also load up the objects we have already inspected and passed judgement upon.
# That information lives here:

lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"
inspect_ngc = pd.read_csv(
    os.path.join(lc_dir, "inspection_ngc.csv"), skipinitialspace=True
)
inspect_ic = pd.read_csv(
    os.path.join(lc_dir, "inspection_ic.csv"), skipinitialspace=True
)

rejected_sources_ngc = inspect_ngc["SOURCEID"][inspect_ngc["exclude?"] == "yes"]
rejected_sources_ic = inspect_ic["SOURCEID"][inspect_ic["exclude?"] == "yes"]

approved_sources_ngc = inspect_ngc["SOURCEID"][inspect_ngc["exclude?"] != "yes"]
approved_sources_ic = inspect_ic["SOURCEID"][inspect_ic["exclude?"] != "yes"]

# and the "script" that created those lightcurves is the following Jupyter notebook:
# /Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/Lightcurves of all matches - first pass for inspection.ipynb


# Basically, we want to identify all the objects that need to be evaluated, and
# generate light curves for them.

names = ["ngc", "ic"]
longnames_dict = {"ngc": "NGC 1333", "ic": "IC 348"}
matches = [ngc_match, ic_match]
inspects = [inspect_ngc, inspect_ic]

photometry_dict = {
    "ngc": photometry.group_wserv_v2(photometry.load_wserv_v2(7)),
    "ic": photometry.group_wserv_v2(photometry.load_wserv_v2(8)),
}

quality_dict = {
    "ngc": quality_classes.load_q(7),
    "ic": quality_classes.load_q(8),
}

lc_fn_dict = {
    "ngc": ngc1333_simple_lc_scatter_brokenaxes,
    "ic": ic348_simple_lc_scatter_brokenaxes,
}

spread_dict = {"ngc": spreadsheet.load_wserv_v2(7), "ic": spreadsheet.load_wserv_v2(8)}

for name, match, inspect in zip(names, matches, inspects):

    spread = spread_dict[name]
    q = quality_dict[name]

    under_consideration = match.all_matches["SpT"] >= 0
    unchecked = ~np.in1d(match.all_matches["SOURCEID"], inspect["SOURCEID"])

    gonna_check = unchecked & under_consideration
    print(
        name,
        "Gonna check",
        np.sum(gonna_check),
        " of ",
        np.sum(under_consideration),
        " under consideration",
    )
    print(
        name, "Already checked", np.sum(~unchecked), "... should equal ", len(inspect)
    )

    # print(match.all_matches[gonna_check]["SOURCEID"])

    gonna_check_sids = match.all_matches[gonna_check]["SOURCEID"]
    gonna_check_table = match.all_matches[gonna_check]

    new_lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3.5"

    # switch to True if you need to make a new table. This overwrites...
    if False:
        with open(os.path.join(new_lc_dir, f"new_inspection_{name}.csv"), "w+") as f:

            f.write("i, Name, SOURCEID, ")

            f.write("J, ")
            f.write("H, ")
            f.write("K, ")

            f.write("Q,  ")
            f.write("S,  ")
            f.write("exclude?,  ")
            f.write("notes,  ")

            f.write("\n")
            for i, sid in enumerate(gonna_check_sids):

                f.write(f"{i}, {gonna_check_table['Name'][i]}, {sid}, ")

                f.write(f"{spread['median']['JAPERMAG3'][sid]}, ")
                f.write(f"{spread['median']['HAPERMAG3'][sid]}, ")
                f.write(f"{spread['median']['KAPERMAG3'][sid]}, ")

                f.write(f"Q{q_string(sid, spread, q):4s},  ")
                f.write(f"{spread['variability']['Stetson_JHK'][sid]},  ")

                f.write("  ,")
                f.write("  ,")

                f.write("\n")

    dat = photometry_dict[name]
    longname = longnames_dict[name]
    lc_fn = lc_fn_dict[name]

    # This section produces the lightcurves.
    for i, sid in enumerate(gonna_check_sids):

        suptitle = ""

        suptitle += (
            "\n"
            + f"{i} {longname} : {gonna_check_table['Name'][i]} \n is {sidsep(sid)}."
        )
        suptitle += "\n" + f" Match: {gonna_check_table['d2d_arcsec'][i]:.2f} arcsec"
        suptitle += "\n" + f" Published spectral type: {gonna_check_table['Adopt'][i]}"
        suptitle += (
            "\n"
            + f" Stetson variability index:  S={spread['variability']['Stetson_JHK'][sid]:.2f}"
        )
        suptitle += "\n" + " "

        fig_lc = lc_fn(dat, sid, cmap="jet")
        #     fig_lc.ax_j.set_title(f"{i} S-{table3['SONYC'][i]}. $S = {spread.wserv7['variability']['Stetson_JHK'].values[idx3[i]]:.2f}$")
        fig_lc.text(
            0.1,
            0.89,
            suptitle,
            horizontalalignment="left",
            transform=fig_lc.transFigure,
            font="monospace",
            fontsize=14,
        )

        data_quality_text = f"Q{q_string(sid, spread, q):4s}  "

        data_quality_text += f"N_J= {int(spread['count']['N_J'][sid]):3d} |"
        data_quality_text += f" {int(spread['count']['N_J_good'][sid]):3d}"
        data_quality_text += f" {int(spread['count']['N_J_info'][sid]):3d}"

        data_quality_text += f"\n       N_H= {int(spread['count']['N_H'][sid]):3d} |"
        data_quality_text += f" {int(spread['count']['N_H_good'][sid]):3d}"
        data_quality_text += f" {int(spread['count']['N_H_info'][sid]):3d}"

        data_quality_text += f"\n       N_K= {int(spread['count']['N_K'][sid]):3d} |"
        data_quality_text += f" {int(spread['count']['N_K_good'][sid]):3d}"
        data_quality_text += f" {int(spread['count']['N_K_info'][sid]):3d}"
        data_quality_text += f"\n       Pstar={spread['median']['PSTAR'][sid]:.2f}"

        #     print(data_quality_text)

        fig_lc.text(
            0.65,
            1.095,
            data_quality_text,
            horizontalalignment="left",
            verticalalignment="top",
            transform=fig_lc.transFigure,
            font="monospace",
            fontsize=14,
        )

        # break
        # if i > 10:
        #     break

        fig_lc.savefig(
            os.path.join(new_lc_dir, f"{name}_{i:03d}_{sid}.png"),
            bbox_inches="tight",
            dpi=100,
        )
        plt.close(fig_lc)
