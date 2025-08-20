"""
In this script I am coding up maps of the regions we have studied. 

Copied from `map_figure_Perseus_regions.py`.

"""

import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wuvars.analysis.variability_selection as sv
from wuvars.analysis.bd_matching_ic348 import lowmass_ic_joint_matches
from wuvars.analysis.bd_matching_ngc1333 import lowmass_ngc_joint_matches
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.spectral_type_to_number import get_SpT_from_num
from wuvars.analysis.variability_selection_curved import (curve_Stetson, sv_hk,
                                                          sv_jh, sv_jhk, sv_jk)
from wuvars.data import photometry, quality_classes, spreadsheet

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

    # ax_ic.plot(
    #     np.degrees(ic_spread["median"]["RA"]),
    #     np.degrees(ic_spread["median"]["DEC"]),
    #     "k,",
    #     alpha=0.25,
    #     rasterized=True,
    # )

    # ax_ic.plot(
    #     np.degrees(lowmass_ic_joint_matches["RA"][~rejected_sources_ic_ind]),
    #     np.degrees(lowmass_ic_joint_matches["DEC"][~rejected_sources_ic_ind]),
    #     "k.",
    #     ms=4,
    #     label="Good match (n=225)",
    # )
    # ax_ic.plot(
    #     np.degrees(lowmass_ic_joint_matches["RA"][rejected_sources_ic_ind]),
    #     np.degrees(lowmass_ic_joint_matches["DEC"][rejected_sources_ic_ind]),
    #     "rx",
    #     ms=6,
    #     mew=1.2,
    #     label="Unusable lightcurve (n=9)",
    # )
    # ax_ic.invert_xaxis()
    # ax_ic.set_xlabel("RA (deg)")
    # ax_ic.set_ylabel("Dec (deg)")

    # ax_ic.set_xlim(
    #     np.degrees(ic_spread["median"]["RA"][ic_q.q2]).max(),
    #     np.degrees(ic_spread["median"]["RA"][ic_q.q2]).min(),
    # )
    # ax_ic.set_ylim(
    #     np.degrees(ic_spread["median"]["DEC"]).min(),
    #     np.degrees(ic_spread["median"]["DEC"]).max(),
    # )

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


if __name__ == "__main__":
    if True:
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
            np.degrees(ngc_match.approved['median_RA']),
            np.degrees(ngc_match.approved['median_DEC']),
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
