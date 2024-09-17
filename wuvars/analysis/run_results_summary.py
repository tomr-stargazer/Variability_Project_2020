"""
This is a script that runs through and summarizes what's up with our data.

"""
import warnings

import matplotlib.pyplot as plt
import numpy as np
import wuvars.analysis.variability_selection as sv
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.spectral_type_to_number import get_SpT_from_num
from wuvars.analysis.variability_selection_curved import (curve_Stetson, sv_hk,
                                                          sv_jh, sv_jhk, sv_jk)
from wuvars.data import photometry, quality_classes, spreadsheet

warnings.filterwarnings("ignore")

# This part pulls up the source information from our position-matching scheme.
ngc_match = match_ngc()
ic_match = match_ic()

names = ["ngc", "ic"]
wserv_dict = {"ngc": 7, "ic": 8}
fullname_dict = {"ngc": "NGC 1333", "ic": "IC 348"}
match_dict = {"ngc": ngc_match, "ic": ic_match}
spread_dict = {"ngc": spreadsheet.load_wserv_v2(7), "ic": spreadsheet.load_wserv_v2(8)}
q_dict = {"ngc": quality_classes.load_q(7), "ic": quality_classes.load_q(8)}

make_figs = True

# Let's do the overview of IC 348 and NGC 1333.
# How many objects were in the original catalog
print("The input catalogs (total):\n")
for name in names:
    print("***")
    print(fullname_dict[name], "\n")
    match = match_dict[name]
    spread = spread_dict[name]

    # Simple stats on the input catalog
    print(f"{len(match.input_catalog)} objects in the input catalog")
    valid = match.input_catalog["SpT"] >= 0.0
    print(f"{np.sum(valid)} with an assigned spectral type of M0 or later")
    print("   and are therefore considered for this study.")
    invalid_nan = np.isnan(match.input_catalog["SpT"])
    invalid_early = match.input_catalog["SpT"] < 0.0
    print(f"({np.sum(invalid_nan)} have no assigned spectral type)")
    print(f"({np.sum(invalid_early)} have an assigned spectral type earlier than M0)")
    assert len(match.input_catalog) == np.sum(valid) + np.sum(invalid_nan) + np.sum(
        invalid_early
    )

    print("")
    print(
        f"{len(match.ML)} objects ({len(match.ML)/np.sum(valid)*100:.1f}%) had a positional match within 0.5''"
    )
    print(f"({len(match.unmatched)} objects failed to positionally match)")

    print("")
    # print(len(match.approved), "Approved")
    print(
        f"{len(match.approved)} objects' light curves ({len(match.approved)/len(match.ML)*100:.1f}%) were approved by manual inspection"
    )
    print(f"({len(match.rejected)} light curves were rejected by manual inspection)")

    # let's divide up by infrared excess
    ir_exc = match.approved["IRexc"] == "yes"

    print(
        f"{np.sum(ir_exc)} of the {len(match.approved)} approved objects have an IR excess"
        f" ({100 * np.sum(ir_exc) / len(match.approved):.1f}%)"
    )

    print("")
    if make_figs:

        # let's put an HR diagram here
        fig_hr, ax_hr = plt.subplots(figsize=(6, 6))

        ax_hr.plot(
            spread["median"]["HMKPNT"],
            spread["median"]["KAPERMAG3"],
            "k,",
            alpha=0.25,
        )
        ax_hr.plot(
            match.approved["median_HMKPNT"],
            match.approved["median_KAPERMAG3"],
            "k.",
            label="without IR exc",
        )
        ax_hr.plot(
            match.approved["median_HMKPNT"][ir_exc],
            match.approved["median_KAPERMAG3"][ir_exc],
            "r.",
            label="with IR exc",
        )

        ax_hr.set_title(f"CMD in {fullname_dict[name]}")
        ax_hr.set_xlabel("$H-K$ color")
        ax_hr.set_ylabel("$K$ mag")
        ax_hr.invert_yaxis()
        ax_hr.set_xlim(0, None)
        ax_hr.legend()

        fig, ax = plt.subplots(figsize=(6, 3))
        ax.hist(match.approved["SpT"], range=[0, 15], bins=np.arange(0, 15, 0.5))

        ax.set_xlim(0, 14)
        ax.set_ylabel("Number of sources")
        ax.set_xlabel("Spectral Type")

        spt_array = np.array([get_SpT_from_num(int(x)) for x in ax.get_xticks()])
        ax.xaxis.set_tick_params(labelbottom=True)
        ax.set_xticklabels(spt_array)

        ax.text(0.025, 0.9, fullname_dict[name], transform=ax.transAxes)
        ax.text(
            0.025,
            0.825,
            str(len(match.approved)) + " approved objects",
            transform=ax.transAxes,
        )

        # ax.axvline(4.7, color="k", ls=":")
        # ax.axvline(6.8, color="k", ls=":")
        ax.axvspan(4.7, 6.8, hatch="xxxx", facecolor="None", ec="k", lw=0.25, alpha=0.2)
        ax.grid(True, axis="x", ls=":")
        plt.show()

    # now... data quality (error) versus spectral type

    if make_figs:
        fig, axes = plt.subplots(figsize=(6, 8), nrows=2)
        ax = axes[0]
        ax.plot(match.approved["SpT"], match.approved["median_KAPERMAG3"], "k.")
        ax.invert_yaxis()

        ax.set_xlim(0, 14)
        ax.set_ylabel("median $K$ mag (UKIRT)")
        ax.set_xlabel("Spectral Type")
        ax.grid(True, axis="x", ls=":")
        ax.axvspan(4.7, 6.8, hatch="xxxx", facecolor="None", ec="k", lw=0.25, alpha=0.1)

        ax1 = axes[1]
        ax1.plot(match.approved["SpT"], match.approved["median_KAPERMAG3ERR"], "k.")
        ax1.set_xlim(0, 14)
        ax1.semilogy()
        ax1.set_xlabel("Spectral Type")
        ax1.set_ylabel("median $K$ photometric uncertainty (mag)")

        ax.set_xticklabels(spt_array)
        ax1.set_xticklabels(spt_array)
        ax1.grid(True, axis="x", ls=":")
        ax1.axvspan(
            4.7, 6.8, hatch="xxxx", facecolor="None", ec="k", lw=0.25, alpha=0.1
        )

        plt.show()

    # print("IC 348:")
    # print("NGC 1333:"
    # print(f"{len(ngc_match.input_catalog)} objects in IC 348;")
    # print(f"{} objects in NGC 1333.")
    # How many objects did we prune down to for our positional matching
    # How many objects were rejected on lightcurve analysis
    # What's the breakdown of quality classes among the "remaining" objects
    # What's the breakdown of spectral types

    print("Variability analysis")
    print(f"{len(match.approved)} approved")
    print(f"{len(match.statistical)} statistical")
    print(f"{len(match.color)} color")

    # borrowed from analysis/prototypes/Prototyping final variable list.ipynb
    # also see
    # wuvars/analysis/prototypes/Prototyping all variables (part 3).ipynb

    """
    okay. At present...
    looks like select_stetson_variables() in prototype2...
    is the function I want!
    Basically, I want to figure out how to re-generate the Stetson vs K mag plots
    with the fun little lines in there showing off "2sigma" and "3sigma" lines.
    This might be more work than I'd hoped for.
    And to print out:
    (a) the breakdown of various quality cuts for the approved sample(s)
    (b) how many objects are automatically detected as variable per the criteria
        I had laid out
    """

    from wuvars.analysis.prototype2_of_variability_full_criteria import (
        select_stetson_variables,
        select_periodic_variables,
        select_periodic_variables_experimental,
    )

    from wuvars.analysis.load_periodics_v4 import select_periodic_variables_v4

    wserv = wserv_dict[name]

    # print(f"{name}: Variability stuff from Stetson stuff: ")
    v2, v1 = select_stetson_variables(wserv)
    # print(np.sum(v2), np.sum(v1), np.sum(v2 & v1))
    # print(len(v2))

    # OKAY. I believe that we "select" the v1 variables that are in the "approved"
    # sample using the following syntax:

    approved_v1 = v1[match.approved["SOURCEID"]]
    approved_v2 = v2[match.approved["SOURCEID"]]
    print(
        f"Number of approved sources that are v1 variables: {np.sum(approved_v1)}/{len(match.approved)}"
    )

    av1 = match.approved

    statistical_v1 = v1[match.statistical["SOURCEID"]]
    statistical_v2 = v2[match.statistical["SOURCEID"]]
    print(
        f"Number of statistical sources that are v1 variables: {np.sum(statistical_v1)}/{len(match.statistical)}"
        f" ({100*np.sum(statistical_v1)/len(match.statistical):.2f}%)"
    )

    # trying the new one!
    v_per = select_periodic_variables_v4(wserv)
    # v_per is

    print(
        f"Number of periodic variables found: {len(v_per)} / {len(match.approved)}"
        f" ({100*len(v_per)/len(match.approved):.2f}%)"
    )
    # print(f"Sum of v_per: {np.sum(v_per)}")

    periodics = np.in1d(match.approved["SOURCEID"], v_per.index)

    if make_figs:
        fig2, axes2 = plt.subplots(figsize=(6, 12*4/3), nrows=4, sharex=True)
        ax2, ax2_2, ax2_3, ax2_4 = axes2
        ax2.plot(
            match.approved["SpT"][periodics],
            match.approved["std_KAPERMAG3"][periodics],
            marker="o",
            ls="None",
            mec="k",
            mfc="None",
            ms=8,
            label="periodic variable",
        )

        ax2.plot(
            match.approved["SpT"][approved_v1],
            match.approved["std_KAPERMAG3"][approved_v1],
            "r.",
            ms=5,
            label="automatic Stetson variable",
        )
        ax2.plot(
            match.approved["SpT"],
            match.approved["std_KAPERMAG3"],
            "k.",
            ms=2,
            label="not automatic Stetson variable",
            zorder=-5,
        )
        ax2.set_xlim(-0.2, 14.2)
        ax2.grid(True, axis="x", ls=":")
        ax2.axvspan(
            4.7, 6.8, hatch="xxxx", facecolor="None", ec="k", lw=0.25, alpha=0.1
        )
        ax2.legend(loc="lower right", fontsize=6)

        ax2_2.plot(
            match.approved["SpT"][periodics],
            match.approved["var_Stetson_JHK"][periodics],
            marker="o",
            ls="None",
            mec="k",
            mfc="None",
            ms=8,
            label="periodic variable",
        )

        ax2_2.plot(
            match.approved["SpT"][approved_v1],
            match.approved["var_Stetson_JHK"][approved_v1],
            "r.",
            ms=5,
            label="automatic Stetson variable",
        )
        ax2_2.plot(
            match.approved["SpT"],
            match.approved["var_Stetson_JHK"],
            "k.",
            ms=2,
            label="not automatic Stetson variable",
            zorder=-5,
        )
        ax2_2.set_xlim(-0.2, 14.2)
        ax2_2.grid(True, axis="x", ls=":")
        ax2_2.axvspan(
            4.7, 6.8, hatch="xxxx", facecolor="None", ec="k", lw=0.25, alpha=0.1
        )
        ax2_2.axhline(1, color="k", lw=0.25, ls="--")

        ax2_2.legend(loc="lower right", fontsize=6)

        ###
        ax2_3.plot(
            match.approved["SpT"][periodics],
            match.approved["var_K_red_chisq"][periodics],
            marker="o",
            ls="None",
            mec="k",
            mfc="None",
            ms=8,
            label="periodic variable",
        )

        ax2_3.plot(
            match.approved["SpT"][approved_v1],
            match.approved["var_K_red_chisq"][approved_v1],
            "r.",
            ms=5,
            label="automatic Stetson variable",
        )
        ax2_3.plot(
            match.approved["SpT"],
            match.approved["var_K_red_chisq"],
            "k.",
            ms=2,
            label="not automatic Stetson variable",
            zorder=-5,
        )
        ax2_3.set_xlim(-0.2, 14.2)
        ax2_3.grid(True, axis="x", ls=":")
        ax2_3.axvspan(
            4.7, 6.8, hatch="xxxx", facecolor="None", ec="k", lw=0.25, alpha=0.1
        )
        ax2_3.legend(loc="upper right", fontsize=6)
        ####

        ax2_4.plot(
            match.approved["SpT"][periodics],
            match.approved["median_KAPERMAG3ERR"][periodics],
            marker="o",
            ls="None",
            mec="k",
            mfc="None",
            ms=8,
            label="periodic variable",
        )

        ax2_4.plot(
            match.approved["SpT"][approved_v1],
            match.approved["median_KAPERMAG3ERR"][approved_v1],
            "r.",
            ms=5,
            label="automatic Stetson variable",
        )
        ax2_4.plot(
            match.approved["SpT"],
            match.approved["median_KAPERMAG3ERR"],
            "k.",
            ms=2,
            label="not automatic Stetson variable",
            zorder=-5,
        )
        ax2_4.set_xlim(-0.2, 14.2)
        ax2_4.grid(True, axis="x", ls=":")
        ax2_4.axvspan(
            4.7, 6.8, hatch="xxxx", facecolor="None", ec="k", lw=0.25, alpha=0.1
        )
        ax2_4.legend(loc="upper left", fontsize=6)
        ax2_4.semilogy()
        ax2_4.set_ylabel("median $K$ photometric uncertainty (mag)")


        ####

        spt_array = np.array([get_SpT_from_num(int(x)) for x in ax2.get_xticks()])
        ax2.xaxis.set_tick_params(labelbottom=True)
        ax2.set_xticklabels(spt_array)
        ax2_2.xaxis.set_tick_params(labelbottom=True)
        ax2_2.set_xticklabels(spt_array)

        ax2.semilogy()
        ax2.set_xlabel("Spectral Type")
        ax2.set_ylabel("rms $K$ variability (mag)")

        ax2_2.semilogy()
        ax2_2.set_ylim(1e-2, None)

        ax2_2.set_xlabel("Spectral Type")
        ax2_2.set_ylabel("$JHK$ Stetson index")

        ax2_3.semilogy()
        # ax2_3.set_ylim(1e-2, None)

        ax2_3.set_xlabel("Spectral Type")
        ax2_3.set_ylabel(r"$\chi^2_{\nu}$ ($K$)")

        fig3, axes3 = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)

        axes3[0].plot(
            match.approved["range_KAPERMAG3"][ir_exc],
            match.approved["SpT"][ir_exc],
            "r.",
        )
        axes3[1].plot(
            match.approved["range_KAPERMAG3"][~ir_exc],
            match.approved["SpT"][~ir_exc],
            "k.",
        )
        axes3[0].invert_yaxis()

        fig4, ax4 = plt.subplots(figsize=(8, 5))

        ax4.plot(
            match.approved["range_KAPERMAG3"][ir_exc],
            match.approved["SpT"][ir_exc],
            "r+",
        )
        ax4.plot(
            match.approved["range_KAPERMAG3"][~ir_exc],
            match.approved["SpT"][~ir_exc],
            "k+",
        )
        ax4.invert_yaxis()
        ax4.set_ylim(13.5, -0.5)
        ax4.set_xlim(0, 1.2)
        ax4.set_title(name)

        plt.show()

    # okay. So, let's re-summarize some things...
    # Variables were selected using 3 criteria:
    # 1. the automated Stetson selection (with a magnitude-dependent Stetson cut)
    # 2. the periodic selection (partially manual interpretation of Lomb-Scargle
    #    periodograms + some followup including detrending of secular variables)
    # 3. the “subjective” selection (partially manual inspection of light curves
    #    + Stetson values)

    print("TBD")
    print("")
# Variability finding
# What counts as a variable star?
# Automatic variability finding: print out how many qualify for this
# How many are "automatically detected" variables?
# - How about "the statistical sample" here?
# Periodic (requires some manual intervention)
# - This cannot be summarized until I've done that
# Visual inspection of

### Classifying variability
# Showing the Q, M plot (for variables only)

# Summarizing "types" of variability
# Don't forget color

# How do these "types" vary with

### Statistical (quantitative) variability trends across spectral type
# Period trends
# Amplitude trends
# (etc)
