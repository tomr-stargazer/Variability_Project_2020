"""
In this script I am coding up maps of the regions we have studied. 

I want to set this up for publication.

"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from wuvars.analysis.bd_matching_ic348 import lowmass_ic_joint_matches
from wuvars.analysis.bd_matching_ngc1333 import lowmass_ngc_joint_matches
from wuvars.data import photometry, quality_classes, spreadsheet

lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"
ngc_spread = spreadsheet.load_wserv_v2(7)
ngc_q = quality_classes.load_q(7)
inspect_ngc = pd.read_csv(
    os.path.join(lc_dir, "inspection_ngc.csv"), skipinitialspace=True
)

rejected_sources_ngc = inspect_ngc["SOURCEID"][inspect_ngc["exclude?"] == "yes"]
# inspect_onc['SOURCEID']

rejected_sources_ngc_ind = np.in1d(
    lowmass_ngc_joint_matches["SOURCEID"], rejected_sources_ngc
)

# Step one: reproduce literally what has been done.


def legacy_NGC1333_map_figure():

    # I am using the following file:
    # ./wuvars/analysis/prototypes/Lightcurves of all matches - first pass for inspection.ipynb

    plt.figure(figsize=(8, 8), dpi=150)

    plt.plot(
        np.degrees(ngc_spread["median"]["RA"]),
        np.degrees(ngc_spread["median"]["DEC"]),
        "k,",
        alpha=0.5,
        rasterized=True,
    )

    plt.plot(
        np.degrees(lowmass_ngc_joint_matches["RA"][~rejected_sources_ngc_ind]),
        np.degrees(lowmass_ngc_joint_matches["DEC"][~rejected_sources_ngc_ind]),
        "k.",
        label="Luhman LMS+BDs with usable lightcurve data (n=103)",
    )
    plt.plot(
        np.degrees(lowmass_ngc_joint_matches["RA"][rejected_sources_ngc_ind]),
        np.degrees(lowmass_ngc_joint_matches["DEC"][rejected_sources_ngc_ind]),
        "rx",
        ms=6,
        mew=1.7,
        label="Luhman LMS+BDs with unusable lightcurve data (n=2)",
    )
    plt.gca().invert_xaxis()
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")

    plt.xlim(
        np.degrees(ngc_spread["median"]["RA"][ngc_q.q2]).max(),
        np.degrees(ngc_spread["median"]["RA"][ngc_q.q2]).min(),
    )
    plt.ylim(
        np.degrees(ngc_spread["median"]["DEC"]).min(),
        np.degrees(ngc_spread["median"]["DEC"]).max(),
    )

    plt.legend(fontsize=14, loc="lower left")
    plt.title("Accepted and rejected source matches in NGC 1333")

    plt.savefig("test_NGC1333_map.png")
    plt.savefig("test_NGC1333_map.pdf")


# === IC 348 ===


lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"
inspect_ic = pd.read_csv(
    os.path.join(lc_dir, "inspection_ic.csv"), skipinitialspace=True
)

ic_spread = spreadsheet.load_wserv_v2(8)
ic_q = quality_classes.load_q(8)

rejected_sources_ic = inspect_ic["SOURCEID"][inspect_ic["exclude?"] == "yes"]
# inspect_onc['SOURCEID']

rejected_sources_ic_ind = np.in1d(
    lowmass_ic_joint_matches["SOURCEID"], rejected_sources_ic
)


def legacy_IC348_map_figure():
    plt.figure(figsize=(8, 8), dpi=150)

    plt.plot(
        np.degrees(ic_spread["median"]["RA"]),
        np.degrees(ic_spread["median"]["DEC"]),
        "k,",
        alpha=0.5,
        rasterized=True,
    )

    plt.plot(
        np.degrees(lowmass_ic_joint_matches["RA"][~rejected_sources_ic_ind]),
        np.degrees(lowmass_ic_joint_matches["DEC"][~rejected_sources_ic_ind]),
        "k.",
        label="Luhman LMS+BDs with usable lightcurve data (n=225)",
    )
    plt.plot(
        np.degrees(lowmass_ic_joint_matches["RA"][rejected_sources_ic_ind]),
        np.degrees(lowmass_ic_joint_matches["DEC"][rejected_sources_ic_ind]),
        "rx",
        ms=6,
        mew=1.7,
        label="Luhman LMS+BDs with unusable lightcurve data (n=9)",
    )
    plt.gca().invert_xaxis()
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")

    plt.xlim(
        np.degrees(ic_spread["median"]["RA"][ic_q.q2]).max(),
        np.degrees(ic_spread["median"]["RA"][ic_q.q2]).min(),
    )
    plt.ylim(
        np.degrees(ic_spread["median"]["DEC"]).min(),
        np.degrees(ic_spread["median"]["DEC"]).max(),
    )

    plt.legend(fontsize=14, loc="lower left")
    plt.title("Accepted and rejected source matches in IC 348")

    plt.savefig("test_IC348_map.png")
    plt.savefig("test_IC348_map.pdf")


def combined_legacy_maps():
    pass


if __name__ == "__main__":
    if False:
        legacy_NGC1333_map_figure()
        legacy_IC348_map_figure()

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
        np.degrees(lowmass_ngc_joint_matches["RA"][~rejected_sources_ngc_ind]),
        np.degrees(lowmass_ngc_joint_matches["DEC"][~rejected_sources_ngc_ind]),
        "k.",
        ms=4,
        label="Good match (n=103)",
    )
    ax_ngc.plot(
        np.degrees(lowmass_ngc_joint_matches["RA"][rejected_sources_ngc_ind]),
        np.degrees(lowmass_ngc_joint_matches["DEC"][rejected_sources_ngc_ind]),
        "rx",
        ms=6,
        mew=1.2,
        label="Unusable lightcurve (n=2)",
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
        np.degrees(lowmass_ic_joint_matches["RA"][~rejected_sources_ic_ind]),
        np.degrees(lowmass_ic_joint_matches["DEC"][~rejected_sources_ic_ind]),
        "k.",
        ms=4,
        label="Good match (n=225)",
    )
    ax_ic.plot(
        np.degrees(lowmass_ic_joint_matches["RA"][rejected_sources_ic_ind]),
        np.degrees(lowmass_ic_joint_matches["DEC"][rejected_sources_ic_ind]),
        "rx",
        ms=6,
        mew=1.2,
        label="Unusable lightcurve (n=9)",
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

    plt.savefig("test_IC348_NGC1333_map.png", bbox_inches='tight')
    plt.savefig("test_IC348_NGC1333_map.pdf", bbox_inches='tight')
