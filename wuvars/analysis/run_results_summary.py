"""
This is a script that runs through and summarizes what's up with our data.

"""
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
qm_dict = {
    "ngc": pd.read_hdf(os.path.join(os.getcwd(), "ngc_QM.h5")).set_index("SOURCEID"),
    "ic": pd.read_hdf(os.path.join(os.getcwd(), "ic_QM.h5")).set_index("SOURCEID"),
}

make_figs = True


if __name__ == "__main__":
    # Let's do the overview of IC 348 and NGC 1333.
    # How many objects were in the original catalog
    print("The input catalogs (total):\n")
    for name in names:
        print("***")
        print(fullname_dict[name], "\n")
        match = match_dict[name]
        spread = spread_dict[name]
        qm = qm_dict[name]

        # Simple stats on the input catalog
        print(f"{len(match.input_catalog)} objects in the input catalog")
        valid = match.input_catalog["SpT"] >= 0.0
        print(f"{np.sum(valid)} with an assigned spectral type of M0 or later")
        print("   and are therefore considered for this study.")
        invalid_nan = np.isnan(match.input_catalog["SpT"])
        invalid_early = match.input_catalog["SpT"] < 0.0
        print(f"({np.sum(invalid_nan)} have no assigned spectral type)")
        print(
            f"({np.sum(invalid_early)} have an assigned spectral type earlier than M0)"
        )
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
        print(
            f"({len(match.rejected)} light curves were rejected by manual inspection)"
        )

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
            ax.hist(
                match.statistical["SpT"],
                edgecolor="k",
                facecolor="None",
                hatch="..",
                range=[0, 15],
                bins=np.arange(0, 15, 0.5),
                histtype="stepfilled",
                label="'statistical' sample\n($K_{\\rm{err}} \\leq 0.05$ mag)\n"
                f"$n$={len(match.statistical)}",
            )

            ax.set_xlim(0, 14)
            ax.set_ylabel("Number of sources")
            ax.set_xlabel("Spectral Type")
            ax.legend()

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
            ax.axvspan(
                4.7, 6.8, hatch="xxxx", facecolor="None", ec="k", lw=0.25, alpha=0.2
            )
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
            ax.axvspan(
                4.7, 6.8, hatch="xxxx", facecolor="None", ec="k", lw=0.25, alpha=0.1
            )

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

        print(f"{fullname_dict[name]} Variability analysis")
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

        # now let's select the subjectives
        from wuvars.analysis.load_subjective_variables import (
            select_subjective_variables,
        )

        v_subj = select_subjective_variables(wserv)

        print(
            f"Number of subjective variables found: {len(v_subj)} / {len(match.approved)}"
            f" ({100*len(v_subj)/len(match.approved):.2f}%)"
        )

        subjectives = np.in1d(match.approved["SOURCEID"], v_subj.index)

        # Final count of variables
        # and nonvariables

        # Total number of objects designated as variable

        vars = approved_v1 | periodics | subjectives
        print(f"Total number of variables: {np.sum(vars)}/{len(vars)}")
        print(f" ({100*np.sum(vars)/len(match.approved):.2f}%)")

        # Total number of objects not designated as variable
        print(f"Total number of not variables: {np.sum(~vars)}/{len(vars)}")
        print(f" ({100*np.sum(~vars)/len(match.approved):.2f}%)")

        # Now I want Total variability rate among “statistical sample”
        # and
        # Periodic variability rate among “statistical sample”

        statisticals = np.in1d(
            match.approved["SOURCEID"], match.statistical["SOURCEID"]
        )

        print(f"{np.sum(vars & statisticals)}")

        # let's make figures.

        if make_figs:
            fig2, axes2 = plt.subplots(figsize=(6, 12 * 4 / 3), nrows=4, sharex=True)
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
                match.approved["SpT"][subjectives],
                match.approved["std_KAPERMAG3"][subjectives],
                "s",
                c="#800020",
                ms=2,
                label="subjective variable",
            )
            ax2.plot(
                match.approved["SpT"],
                match.approved["std_KAPERMAG3"],
                "k.",
                ms=2,
                label="not identified as variable",
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
                match.approved["SpT"][subjectives],
                match.approved["var_Stetson_JHK"][subjectives],
                "s",
                c="#800020",
                ms=2,
                label="subjective variable",
            )

            ax2_2.plot(
                match.approved["SpT"],
                match.approved["var_Stetson_JHK"],
                "k.",
                ms=2,
                label="not identified as variable",
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
                match.approved["SpT"][subjectives],
                match.approved["var_K_red_chisq"][subjectives],
                "s",
                c="#800020",
                ms=2,
                label="subjective variable",
            )

            ax2_3.plot(
                match.approved["SpT"],
                match.approved["var_K_red_chisq"],
                "k.",
                ms=2,
                label="not identified as variable",
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
                match.approved["SpT"][subjectives],
                match.approved["median_KAPERMAG3ERR"][subjectives],
                "s",
                c="#800020",
                ms=2,
                label="subjective variable",
            )
            ax2_4.plot(
                match.approved["SpT"],
                match.approved["median_KAPERMAG3ERR"],
                "k.",
                ms=2,
                label="not identified as variable",
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

            plot_SpT_vs_rms = False
            if plot_SpT_vs_rms:
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

            # Counting number of good lightcurves to analyze...
            # in each band.

            fig_bands, axes_bands = plt.subplots(nrows=3, ncols=1, sharex=True)

            q = q_dict[name]

            for ax, band in zip(axes_bands, ["J", "H", "K"]):

                approved_q1_b = q[f"q1_{band.lower()}"][match.approved["SOURCEID"]]

                ax.plot(
                    match.approved["SpT"][approved_q1_b],
                    match.approved[f"median_{band}APERMAG3"][approved_q1_b],
                    "k.",
                )
                ax.invert_yaxis()
                ax.text(
                    10,
                    12.5,
                    f"n (>M6) ={np.sum(approved_q1_b & (match.approved['SpT'] > 6))}",
                )
                ax.set_ylabel(band)
                ax.set_xlabel("SpT")
                print("XXX " + band)

                # ok - from this I learned that K band is best for the late types.
                # (unsurprising) but it's good to confirm.

                # Therefore I'll stick to K (and use Q1_K) for the Q, M plots.

                # Building that in soon...

            ####################
            ####################
            ####################
            ####################
            # Here's my workspace.
            # I have been making Q, M plots (sometimes with three-color info)
            # by sticking to the v2 variables.
            # I'd like to expand this to ~all variables~ that are at least q1_k.
            # In other words: their K band is 'good' in a meaningful sense and
            # they are known as variables, either:
            # - stetson/automatic
            # - periodic
            # - subjective

            # Because I'm gonna show their Q_K and M_K stats, it doesn't much
            # matter whether their variability was ~found~ at K band. Just that
            # they ~are variable~ and ~have valid data at K~.

            # So. How do I `select all variables`?
            # combine the selections, right?

            variable = periodics | subjectives | approved_v1
            q1_k = q["q1_k"]
            q1_h = q["q1_h"]            
            q1_j = q["q1_j"]                        
            approved_q1_k = q["q1_k"][match.approved["SOURCEID"]]
            approved_q1_h = q["q1_h"][match.approved["SOURCEID"]]            
            approved_q1_j = q["q1_j"][match.approved["SOURCEID"]]            
            ir_exc_sid = match.approved["SOURCEID"][ir_exc]

            new_fig_j, new_ax_j = plt.subplots(figsize=(6,6))

            new_ax_j.scatter(
                qm["Q_J"][approved_q1_j & variable],
                qm["M_J"][approved_q1_j & variable],
                c='k',
                s=10 + match.approved["range_JAPERMAG3"][approved_q1_j & variable] * 100,
                ec="k",
                lw=0.5,
            )

            new_ax_j.scatter(
                qm["Q_J"][approved_q1_j & variable & ir_exc.data],
                qm["M_J"][approved_q1_j & variable & ir_exc.data],
                c='r',
                s=10 + match.approved["range_JAPERMAG3"][ir_exc.data & approved_q1_j & variable] * 100,
                ec="k",
                lw=0.5,
            )


            new_ax_j.scatter(
                qm["Q_J"][v_per.index][q1_j],
                qm["M_J"][v_per.index][q1_j],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_JAPERMAG3"][approved_q1_j & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )

            new_ax_j.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
            new_ax_j.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)

            new_ax_j.axvline(0.3, 0, 1, color="k", ls=":", lw=0.5)
            new_ax_j.axvline(0.7, 0, 1, color="k", ls=":", lw=0.5)

            new_ax_j.set_xlabel("$Q_J$ score")
            new_ax_j.set_ylabel("$M_J$ score")
            new_ax_j.set_title(
                f"$M$ vs. $Q$ plot for q1_j variables in {fullname_dict[name]}"
            )            

            new_ax_j.set_xlim(-0.3, 0.98)
            new_ax_j.set_ylim(1.275, -1.275)



            new_fig_h, new_ax_h = plt.subplots(figsize=(6,6))

            # for condition, color in zip([ir_exc[approved_q1_k], ~ir_exc[approved_q1_k]], ['k', 'r']):
            #     new_ax_h.scatter(
            #         qm["Q_K"][q1_k][condition],
            #         qm["M_K"][q1_k][condition],
            #         c=color,
            #         s=10 + match.approved["range_KAPERMAG3"][approved_q1_k][condition] * 100,
            #         ec="k",
            #         lw=0.5,
            #     )


            new_ax_h.scatter(
                qm["Q_H"][approved_q1_h & variable],
                qm["M_H"][approved_q1_h & variable],
                c='k',
                s=10 + match.approved["range_HAPERMAG3"][approved_q1_h & variable] * 100,
                ec="k",
                lw=0.5,
            )

            new_ax_h.scatter(
                qm["Q_H"][approved_q1_h & variable & ir_exc.data],
                qm["M_H"][approved_q1_h & variable & ir_exc.data],
                c='r',
                s=10 + match.approved["range_HAPERMAG3"][ir_exc.data & approved_q1_h & variable] * 100,
                ec="k",
                lw=0.5,
            )


            new_ax_h.scatter(
                qm["Q_H"][v_per.index][q1_h],
                qm["M_H"][v_per.index][q1_h],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_HAPERMAG3"][approved_q1_h & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )

            new_ax_h.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
            new_ax_h.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)

            new_ax_h.axvline(0.3, 0, 1, color="k", ls=":", lw=0.5)
            new_ax_h.axvline(0.7, 0, 1, color="k", ls=":", lw=0.5)

            new_ax_h.set_xlabel("$Q_H$ score")
            new_ax_h.set_ylabel("$M_H$ score")
            new_ax_h.set_title(
                f"$M$ vs. $Q$ plot for q1_h variables in {fullname_dict[name]}"
            )            

            new_ax_h.set_xlim(-0.3, 0.98)
            new_ax_h.set_ylim(1.275, -1.275)

            new_fig, new_ax = plt.subplots(figsize=(6,6))

            # for condition, color in zip([ir_exc[approved_q1_k], ~ir_exc[approved_q1_k]], ['k', 'r']):
            #     new_ax.scatter(
            #         qm["Q_K"][q1_k][condition],
            #         qm["M_K"][q1_k][condition],
            #         c=color,
            #         s=10 + match.approved["range_KAPERMAG3"][approved_q1_k][condition] * 100,
            #         ec="k",
            #         lw=0.5,
            #     )


            # TODO: what's a criterion that's Q1_K AND VARIABLE? we need that. Don't include
            # nonvariables in these plots ...
            new_ax.scatter(
                qm["Q_K"][approved_q1_k & variable],
                qm["M_K"][approved_q1_k & variable],
                c='k',
                s=10 + match.approved["range_KAPERMAG3"][approved_q1_k & variable] * 100,
                ec="k",
                lw=0.5,
            )

            new_ax.scatter(
                qm["Q_K"][approved_q1_k & variable & ir_exc.data],
                qm["M_K"][approved_q1_k & variable & ir_exc.data],
                c='r',
                s=10 + match.approved["range_KAPERMAG3"][ir_exc.data & approved_q1_k & variable] * 100,
                ec="k",
                lw=0.5,
            )


            new_ax.scatter(
                qm["Q_K"][v_per.index][q1_k],
                qm["M_K"][v_per.index][q1_k],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_KAPERMAG3"][approved_q1_k & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )

            new_ax.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
            new_ax.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)

            new_ax.axvline(0.3, 0, 1, color="k", ls=":", lw=0.5)
            new_ax.axvline(0.7, 0, 1, color="k", ls=":", lw=0.5)

            new_ax.set_xlabel("$Q_K$ score")
            new_ax.set_ylabel("$M_K$ score")
            new_ax.set_title(
                f"$M$ vs. $Q$ plot for q1_k variables in {fullname_dict[name]}"
            )            

            new_ax.set_xlim(-0.3, 0.98)
            new_ax.set_ylim(1.275, -1.275)


            # Okay. Now we are going to make six additional plots:
            # M vs SpT (J, H, K) 
            # Q vs SpT (J, H, K)
            # these plots will use the ssame coding for ir exc, periodicity, and amplitude as the above.

            fig_Q_vs_M, axes_Q_vs_M = plt.subplots(nrows=3, sharex=True, sharey=True, figsize=(5,10))

            for ax, band in zip(axes_Q_vs_M, ['J', 'H', 'K']):

                q1_b = q[f"q1_{band.lower()}"]
                approved_q1_b = q1_b[match.approved["SOURCEID"]]

                # draw all the points
                ax.scatter(
                    qm[f"Q_{band}"][approved_q1_b & variable],
                    qm[f"M_{band}"][approved_q1_b & variable],
                    c='k',
                    s=10 + match.approved["range_KAPERMAG3"][approved_q1_b & variable] * 100,
                    ec="k",
                    lw=0.5,
                )                

                # draw just the red ones
                ax.scatter(
                    qm[f"Q_{band}"][approved_q1_b & variable & ir_exc.data],
                    qm[f"M_{band}"][approved_q1_b & variable & ir_exc.data],
                    c='r',
                    s=10 + match.approved[f"range_{band}APERMAG3"][ir_exc.data & approved_q1_b & variable] * 100,
                    ec="k",
                    lw=0.5,
                )                

                # draw the little circles
                ax.scatter(
                    qm[f"Q_{band}"][v_per.index][q1_b],
                    qm[f"M_{band}"][v_per.index][q1_b],
                    ec="k",
                    lw=0.5,
                    c="None",
                    s=80 + match.approved[f"range_{band}APERMAG3"][approved_q1_b & periodics] * 100,
                    zorder=10,
                    alpha=0.5,
                )
                # set the axis labels etc
                ax.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
                ax.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)
                ax.axvline(0.3, 0, 1, color="k", ls=":", lw=0.5)
                ax.axvline(0.7, 0, 1, color="k", ls=":", lw=0.5)


                ax.set_xlabel(f"$Q_{band}$ score")
                ax.set_ylabel(f"$M_{band}$ score")
                ax.set_title(
                    f"$M$ vs. $Q$ plot for q1_{band.lower()} variables in {fullname_dict[name]}",
                    fontsize=8
                )            

                ax.set_xlim(-0.3, 0.98)
                ax.set_ylim(1.275, -1.275)


            fig_M_vs_SpT, axes_M_vs_SpT = plt.subplots(nrows=3, sharex=True, sharey=True, figsize=(5,10))

            for ax, band in zip(axes_M_vs_SpT, ['J', 'H', 'K']):

                q1_b = q[f"q1_{band.lower()}"]
                approved_q1_b = q1_b[match.approved["SOURCEID"]]

                # draw all the points
                ax.scatter(
                    match.approved["SpT"][approved_q1_b & variable],
                    qm[f"M_{band}"][approved_q1_b & variable],
                    c='k',
                    s=10 + match.approved["range_KAPERMAG3"][approved_q1_b & variable] * 100,
                    ec="k",
                    lw=0.5,
                )                

                # draw just the red ones
                ax.scatter(
                    match.approved["SpT"][approved_q1_b & variable & ir_exc.data],
                    qm[f"M_{band}"][approved_q1_b & variable & ir_exc.data],
                    c='r',
                    s=10 + match.approved[f"range_{band}APERMAG3"][ir_exc.data & approved_q1_b & variable] * 100,
                    ec="k",
                    lw=0.5,
                )                

                # draw the little circles
                ax.scatter(
                    match.approved["SpT"][approved_q1_b & periodics],
                    qm[f"M_{band}"][v_per.index][q1_b],
                    ec="k",
                    lw=0.5,
                    c="None",
                    s=80 + match.approved[f"range_{band}APERMAG3"][approved_q1_b & periodics] * 100,
                    zorder=10,
                    alpha=0.5,
                )
                # set the axis labels etc
                ax.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
                ax.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)
                # ax.axvline(0.3, 0, 1, color="k", ls=":", lw=0.5)
                # ax.axvline(0.7, 0, 1, color="k", ls=":", lw=0.5)

                if band == 'K':
                    ax.set_xlabel(f"SpT")

                ax.set_ylabel(f"$M_{band}$ score")
                ax.set_title(
                    f"$M$ vs. SpT plot for q1_{band.lower()} variables in {fullname_dict[name]}",
                    fontsize=8
                )            

                ax.set_xlim(-0.1, 10.1)
                ax.invert_yaxis()

                spt_array = np.array([get_SpT_from_num(int(x)) for x in ax.get_xticks()])
                ax.xaxis.set_tick_params(labelbottom=True)
                ax.set_xticklabels(spt_array)

                # ax.set_ylim(1.275, -1.275)

            fig_Q_vs_SpT, axes_Q_vs_SpT = plt.subplots(nrows=3, sharex=True, sharey=True, figsize=(5,10))

            for ax, band in zip(axes_Q_vs_SpT, ['J', 'H', 'K']):

                q1_b = q[f"q1_{band.lower()}"]
                approved_q1_b = q1_b[match.approved["SOURCEID"]]

                # draw all the points
                ax.scatter(
                    match.approved["SpT"][approved_q1_b & variable],
                    qm[f"Q_{band}"][approved_q1_b & variable],
                    c='k',
                    s=10 + match.approved["range_KAPERMAG3"][approved_q1_b & variable] * 100,
                    ec="k",
                    lw=0.5,
                )                

                # draw just the red ones
                ax.scatter(
                    match.approved["SpT"][approved_q1_b & variable & ir_exc.data],
                    qm[f"Q_{band}"][approved_q1_b & variable & ir_exc.data],
                    c='r',
                    s=10 + match.approved[f"range_{band}APERMAG3"][ir_exc.data & approved_q1_b & variable] * 100,
                    ec="k",
                    lw=0.5,
                )                

                # draw the little circles
                ax.scatter(
                    match.approved["SpT"][approved_q1_b & periodics],
                    qm[f"Q_{band}"][v_per.index][q1_b],
                    ec="k",
                    lw=0.5,
                    c="None",
                    s=80 + match.approved[f"range_{band}APERMAG3"][approved_q1_b & periodics] * 100,
                    zorder=10,
                    alpha=0.5,
                )
                # set the axis labels etc
                # ax.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
                # ax.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)
                ax.axhline(0.3, 0, 1, color="k", ls=":", lw=0.5)
                ax.axhline(0.7, 0, 1, color="k", ls=":", lw=0.5)

                if band == 'K':
                    ax.set_xlabel(f"SpT")

                ax.set_ylabel(f"$Q_{band}$ score")
                ax.set_title(
                    f"$Q$ vs. SpT plot for q1_{band.lower()} variables in {fullname_dict[name]}",
                    fontsize=8
                )            

                ax.set_xlim(-0.1, 10.1)
                ax.set_ylim(-0.3, 0.98)
                # ax.invert_yaxis()

                spt_array = np.array([get_SpT_from_num(int(x)) for x in ax.get_xticks()])
                ax.xaxis.set_tick_params(labelbottom=True)
                ax.set_xticklabels(spt_array)


            # Exploring Q, M

            vanilla_v_per = v_per[v_per["Method"] != "poly4"]
            poly_v_per = v_per[v_per["Method"] == "poly4"]

            # first test - can I plot ~anything~
            fig_qm1, ax_qm1 = plt.subplots(figsize=(6, 6))
            ax_qm1.plot(qm["Q_K"][v2], qm["M_K"][v2], "k.")

            ax_qm1.plot(
                qm["Q_K"][poly_v_per.index],
                qm["M_K"][poly_v_per.index],
                marker="s",
                ls="None",
                mec="g",
                mfc="None",
                ms=6,
            )

            ax_qm1.plot(
                qm["Q_K"][vanilla_v_per.index],
                qm["M_K"][vanilla_v_per.index],
                marker="o",
                ls="None",
                mec="k",
                mfc="None",
                ms=8,
            )

            # ax_qm1.set_xlim(-0.5, 1.5)

            ax_qm1.set_title("Prototyping something here.")

            # Okay. Some observations.
            # I expected the periodic circles to be HIGHLY concentrated around 0.
            # Not.. spread evenly throughout.
            # Makes me worried that I am either
            # (a) making this plot wrog (indexing errors etc), or
            # (b) I have computed Q wrong.

            # So - before we proceed in doing anything with Q, we NEED to look at Q scores
            # from various light curves.

            # OK! Did that. I am no longer worried about the Q scores! or my indexing.
            # The culprit was the fact that my period recovery with detrending was ~quite good actually~.

            # build a new one...
            fig_qm1b, ax_qm1b = plt.subplots(figsize=(8, 6))

            sc = ax_qm1b.scatter(
                qm["Q_H"][v2],
                qm["M_H"][v2],
                c=(qm["M_J"][v2] - qm["M_K"][v2]),
                s=10 + match.approved["range_HAPERMAG3"][approved_v2] * 100,
                cmap="RdBu_r",
                vmin=-1,
                vmax=1,
                ec="k",
                lw=0.5,
            )

            ax_qm1b.scatter(
                qm["Q_H"][v_per.index][v2],
                qm["M_H"][v_per.index][v2],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_HAPERMAG3"][approved_v2 & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )

            ax_qm1b.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
            ax_qm1b.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)

            ax_qm1b.axvline(0.3, 0, 1, color="k", ls=":", lw=0.5)
            ax_qm1b.axvline(0.7, 0, 1, color="k", ls=":", lw=0.5)

            plt.colorbar(mappable=sc, label="$M_J - M_K$", ax=ax_qm1b)

            ax_qm1b.set_xlabel("$Q_H$ score")
            ax_qm1b.set_ylabel("$M_H$ score")
            ax_qm1b.set_title(
                f"$M$ vs. $Q$ plot for v2 variables in {fullname_dict[name]}"
            )

            ax_qm1b.set_xlim(-0.4, 0.9)
            ax_qm1b.invert_yaxis()

            # NOW: do the Q, M, SpT one. Color is black to white.
            fig_qm1c, ax_qm1c = plt.subplots(figsize=(8, 6))

            sc = ax_qm1c.scatter(
                qm["Q_H"][v2],
                qm["M_H"][v2],
                c=match.approved["SpT"][approved_v2],
                s=10 + match.approved["range_HAPERMAG3"][approved_v2] * 100,
                cmap="viridis",
                vmin=0,
                vmax=10,
                ec="k",
                lw=0.5,
            )

            ax_qm1c.scatter(
                qm["Q_H"][v_per.index][v2],
                qm["M_H"][v_per.index][v2],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_HAPERMAG3"][approved_v2 & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )

            ax_qm1c.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
            ax_qm1c.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)

            ax_qm1c.axvline(0.3, 0, 1, color="k", ls=":", lw=0.5)
            ax_qm1c.axvline(0.7, 0, 1, color="k", ls=":", lw=0.5)

            plt.colorbar(mappable=sc, label="SpT", ax=ax_qm1c)

            ax_qm1c.set_xlabel("$Q_H$ score")
            ax_qm1c.set_ylabel("$M_H$ score")
            ax_qm1c.set_title(
                f"$M$ vs. $Q$ plot for v2 variables in {fullname_dict[name]}"
            )

            ax_qm1c.set_xlim(-0.4, 0.9)
            ax_qm1c.invert_yaxis()

            # second test -  can I plot stuff from QM versus stuff from approved (e.g. SpT)
            fig_qm2, ax_qm2 = plt.subplots(figsize=(6, 4))

            sc = ax_qm2.scatter(
                match.approved["SpT"][approved_v2],
                qm["M_H"][v2],
                c=(qm["M_J"][v2] - qm["M_K"][v2]),
                s=10 + match.approved["range_HAPERMAG3"][approved_v2] * 100,
                cmap="RdBu_r",
                vmin=-1,
                vmax=1,
                ec="k",
                lw=0.5,
            )
            plt.colorbar(mappable=sc, label="$M_J - M_K$", ax=ax_qm2)

            ax_qm2.scatter(
                match.approved["SpT"][approved_v2 & periodics],
                qm["M_H"][v_per.index][v2],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_HAPERMAG3"][approved_v2 & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )
            ax_qm2.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
            ax_qm2.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)

            ax_qm2.set_xlabel("SpT")
            ax_qm2.set_ylabel("$M_H$ score")
            ax_qm2.set_title(
                f"$M$ vs. SpT plot for v2 variables in {fullname_dict[name]}"
            )

            # ax_qm2.plot(match.approved["SpT"][approved_v2], qm['M_H'][v2], 'k.')

            # # now overlay periodics on this
            # ax_qm2.plot(
            #     match.approved["SpT"][approved_v2 & periodics], qm['M_H'][v_per.index][v2],
            #     marker="o",
            #     ls="None",
            #     mec="k",
            #     mfc="None",
            #     ms=8,)

            ax_qm2.set_xlim(-0.2, 10.5)
            ax_qm2.invert_yaxis()

            qm2_spt_array = np.array(
                [get_SpT_from_num(int(x)) for x in ax_qm2.get_xticks()]
            )
            ax_qm2.xaxis.set_tick_params(labelbottom=True)
            ax_qm2.set_xticklabels(qm2_spt_array)

            fig_qm2b, ax_qm2b = plt.subplots(figsize=(6, 4))

            sc = ax_qm2b.scatter(
                match.approved["SpT"][approved_v2],
                qm["Q_H"][v2],
                c=(qm["M_J"][v2] - qm["M_K"][v2]),
                s=10 + match.approved["range_HAPERMAG3"][approved_v2] * 100,
                cmap="RdBu_r",
                vmin=-1,
                vmax=1,
                ec="k",
                lw=0.5,
            )

            ax_qm2b.scatter(
                match.approved["SpT"][approved_v2 & periodics],
                qm["Q_H"][v_per.index][v2],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_HAPERMAG3"][approved_v2 & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )

            plt.colorbar(mappable=sc, label="$M_J - M_K$", ax=ax_qm2b)

            ax_qm2b.axhline(0.3, 0, 1, color="k", ls=":", lw=0.5)
            ax_qm2b.axhline(0.7, 0, 1, color="k", ls=":", lw=0.5)

            ax_qm2b.set_xlabel("SpT")
            ax_qm2b.set_ylabel("$Q_H$ score")
            ax_qm2b.set_title(
                f"$Q$ vs. SpT plot for v2 variables in {fullname_dict[name]}"
            )

            ax_qm2b.set_xlim(-0.2, 10.5)
            ax_qm2b.set_ylim(-0.3, 0.9)

            qm2b_spt_array = np.array(
                [get_SpT_from_num(int(x)) for x in ax_qm2b.get_xticks()]
            )
            ax_qm2b.xaxis.set_tick_params(labelbottom=True)
            ax_qm2b.set_xticklabels(qm2b_spt_array)

            # third test -  can I plot stuff from QM versus stuff from period
            fig_qm3, ax_qm3 = plt.subplots(figsize=(6, 6))

            plt.show()

            # NOW with ir-excess encoded
            fig_irexc_qm1b, ax_irexc_qm1b = plt.subplots(figsize=(8, 6))

            sc = ax_irexc_qm1b.scatter(
                qm["Q_H"][v2],
                qm["M_H"][v2],
                c=ir_exc[approved_v2],
                s=10 + match.approved["range_HAPERMAG3"][approved_v2] * 100,
                cmap="hot",
                vmin=0,
                vmax=3,
                ec="k",
                lw=0.5,
            )

            ax_irexc_qm1b.scatter(
                qm["Q_H"][v_per.index][v2],
                qm["M_H"][v_per.index][v2],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_HAPERMAG3"][approved_v2 & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )

            ax_irexc_qm1b.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
            ax_irexc_qm1b.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)

            ax_irexc_qm1b.axvline(0.3, 0, 1, color="k", ls=":", lw=0.5)
            ax_irexc_qm1b.axvline(0.7, 0, 1, color="k", ls=":", lw=0.5)

            plt.colorbar(mappable=sc, label="IR_exc (1=yes)", ax=ax_irexc_qm1b)

            ax_irexc_qm1b.set_xlabel("$Q_H$ score")
            ax_irexc_qm1b.set_ylabel("$M_H$ score")
            ax_irexc_qm1b.set_title(
                f"$M$ vs. $Q$ plot for v2 variables in {fullname_dict[name]}"
            )

            ax_irexc_qm1b.set_xlim(-0.4, 0.9)
            ax_irexc_qm1b.invert_yaxis()

            # second test -  can I plot stuff from QM versus stuff from approved (e.g. SpT)
            fig_irexc_qm2, ax_irexc_qm2 = plt.subplots(figsize=(6, 4))

            sc = ax_irexc_qm2.scatter(
                match.approved["SpT"][approved_v2],
                qm["M_H"][v2],
                c=ir_exc[approved_v2],
                s=10 + match.approved["range_HAPERMAG3"][approved_v2] * 100,
                cmap="hot",
                vmin=0,
                vmax=3,
                ec="k",
                lw=0.5,
            )
            plt.colorbar(mappable=sc, label="IR_exc (1=yes)", ax=ax_irexc_qm2)

            ax_irexc_qm2.scatter(
                match.approved["SpT"][approved_v2 & periodics],
                qm["M_H"][v_per.index][v2],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_HAPERMAG3"][approved_v2 & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )
            ax_irexc_qm2.axhline(-0.25, 0, 1, color="k", ls=":", lw=0.5)
            ax_irexc_qm2.axhline(0.25, 0, 1, color="k", ls=":", lw=0.5)

            ax_irexc_qm2.set_xlabel("SpT")
            ax_irexc_qm2.set_ylabel("$M_H$ score")
            ax_irexc_qm2.set_title(
                f"$M$ vs. SpT plot for v2 variables in {fullname_dict[name]}"
            )

            # ax_irexc_qm2.plot(match.approved["SpT"][approved_v2], qm['M_H'][v2], 'k.')

            # # now overlay periodics on this
            # ax_irexc_qm2.plot(
            #     match.approved["SpT"][approved_v2 & periodics], qm['M_H'][v_per.index][v2],
            #     marker="o",
            #     ls="None",
            #     mec="k",
            #     mfc="None",
            #     ms=8,)

            ax_irexc_qm2.set_xlim(-0.2, 10.5)
            ax_irexc_qm2.invert_yaxis()

            qm2_spt_array = np.array(
                [get_SpT_from_num(int(x)) for x in ax_irexc_qm2.get_xticks()]
            )
            ax_irexc_qm2.xaxis.set_tick_params(labelbottom=True)
            ax_irexc_qm2.set_xticklabels(qm2_spt_array)

            fig_irexc_qm2b, ax_irexc_qm2b = plt.subplots(figsize=(6, 4))

            sc = ax_irexc_qm2b.scatter(
                match.approved["SpT"][approved_v2],
                qm["Q_H"][v2],
                c=ir_exc[approved_v2],
                s=10 + match.approved["range_HAPERMAG3"][approved_v2] * 100,
                cmap="hot",
                vmin=0,
                vmax=3,
                ec="k",
                lw=0.5,
            )

            ax_irexc_qm2b.scatter(
                match.approved["SpT"][approved_v2 & periodics],
                qm["Q_H"][v_per.index][v2],
                ec="k",
                lw=0.5,
                c="None",
                s=80 + match.approved["range_HAPERMAG3"][approved_v2 & periodics] * 100,
                zorder=10,
                alpha=0.5,
            )

            plt.colorbar(mappable=sc, label="IR_exc (1=yes)", ax=ax_irexc_qm2b)

            ax_irexc_qm2b.axhline(0.3, 0, 1, color="k", ls=":", lw=0.5)
            ax_irexc_qm2b.axhline(0.7, 0, 1, color="k", ls=":", lw=0.5)

            ax_irexc_qm2b.set_xlabel("SpT")
            ax_irexc_qm2b.set_ylabel("$Q_H$ score")
            ax_irexc_qm2b.set_title(
                f"$Q$ vs. SpT plot for v2 variables in {fullname_dict[name]}"
            )

            ax_irexc_qm2b.set_xlim(-0.2, 10.5)
            ax_irexc_qm2b.set_ylim(-0.3, 0.9)

            qm2b_spt_array = np.array(
                [get_SpT_from_num(int(x)) for x in ax_irexc_qm2b.get_xticks()]
            )
            ax_irexc_qm2b.xaxis.set_tick_params(labelbottom=True)
            ax_irexc_qm2b.set_xticklabels(qm2b_spt_array)

            # third test -  can I plot stuff from QM versus stuff from period
            fig_irexc_qm3, ax_irexc_qm3 = plt.subplots(figsize=(6, 6))

            plt.show()

        # okay. So, let's re-summarize some things...
        # Variables were selected using 3 criteria:
        # 1. the automated Stetson selection (with a magnitude-dependent Stetson cut)
        # 2. the periodic selection (partially manual interpretation of Lomb-Scargle
        #    periodograms + some followup including detrending of secular variables)
        # 3. the “subjective” selection (partially manual inspection of light curves
        #    + Stetson values)

        # print("TBD")
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
