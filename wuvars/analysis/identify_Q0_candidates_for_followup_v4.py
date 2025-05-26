"""
We are identifying Q0 candidates for followup.

"""

import os
import pdb
import warnings

import astropy.table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wuvars.analysis.run_results_summary as rs
import wuvars.analysis.variability_selection as sv
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.q_string import q_string
from wuvars.analysis.spectral_type_to_number import get_SpT_from_num
from wuvars.analysis.variability_selection_curved import (curve_Stetson, sv_hk,
                                                          sv_jh, sv_jhk, sv_jk)
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.plotting.lightcurve import (ic348_simple_lc_brokenaxes,
                                        ic348_simple_lc_scatter_brokenaxes,
                                        ngc1333_simple_lc_scatter_brokenaxes,
                                        simple_lc, simple_lc_brokenaxes,
                                        simple_lc_scatter_brokenaxes)

warnings.filterwarnings("ignore")


# This part pulls up the source information from our position-matching scheme.
rs.ngc_match
rs.ic_match

rs.names
rs.wserv_dict
rs.fullname_dict
rs.match_dict
rs.spread_dict
rs.q_dict

lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/Results/Q0_followup_v4"
non_lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/Results/nonvariables_v4"


def generate_list_for_followup():
    # Let's do the overview of IC 348 and NGC 1333.
    # How many objects were in the original catalog
    print("The input catalogs (total):\n")
    for name in rs.names:
        print("***")
        print(rs.fullname_dict[name], "\n")
        match = rs.match_dict[name]
        spread = rs.spread_dict[name]

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

        wserv = rs.wserv_dict[name]

        # print(f"{name}: Variability stuff from Stetson stuff: ")
        v2, v1, v0 = select_stetson_variables(wserv, return_v0=True)
        # print(np.sum(v2), np.sum(v1), np.sum(v2 & v1))
        # print(len(v2))

        # OKAY. I believe that we "select" the v1 variables that are in the "approved"
        # sample using the following syntax:

        approved_v0 = v0[match.approved["SOURCEID"]]
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

        # now let's

        print("")

        print(
            "This many objects have been identified as BOTH periodic AND an automatic variable:"
        )

        both_autov_and_per = periodics & approved_v1

        print(f"{np.sum(both_autov_and_per)}")

        print("")

        print(
            "This many objects have been identified as NEITHER periodic NOR an automatic variable:"
        )

        neither_autov_nor_per = ~periodics & ~approved_v1

        print(f"{np.sum(neither_autov_nor_per)}")
        print(
            f" (min: {len(match.approved) - np.sum(periodics) - np.sum(approved_v1)})"
        )

        print(
            "Of those, we are going to identify the following number for V0 followup:"
        )

        v0_cand = neither_autov_nor_per & approved_v0
        nonvar = neither_autov_nor_per & ~approved_v0

        print(f"{np.sum(v0_cand)}")
        print(
            f"Also there are {np.sum(nonvar)} that are neither variable nor a v0 candidate."
        )
        # print(match.approved["SOURCEID"][v0_cand])

        # okay!! let's export that list (with Stetson values ready to go)

        try:
            match.approved[v0_cand.values].write(
                os.path.join(
                    lc_dir, name, f"Q0_objects_for_followup_{name}_2024_Oct_13.h5"
                )
            )
            match.approved[v0_cand.values].write(
                os.path.join(
                    lc_dir, name, f"Q0_objects_for_followup_{name}_2024_Oct_13.csv"
                )
            )
        except OSError as e:
            print(e)
            pass

        try:
            match.approved[nonvar.values].write(
                os.path.join(
                    non_lc_dir,
                    name,
                    f"nonvariable_objects_for_followup_{name}_2024_Oct_13.h5",
                )
            )
            match.approved[nonvar.values].write(
                os.path.join(
                    non_lc_dir,
                    name,
                    f"nonvariable_objects_for_followup_{name}_2024_Oct_13.csv",
                )
            )
        except OSError as e:
            print(e)
            pass

def generate_lightcurves_for_inspection():
    dat_dict = {
        "ngc": photometry.group_wserv_v2(
            photometry.load_wserv_v2(7, suffix="_outlier_cleaned_152")
        ),
        "ic": photometry.group_wserv_v2(photometry.load_wserv_v2(8)),
    }
    plot_dict = {
        "ngc": ngc1333_simple_lc_scatter_brokenaxes,
        "ic": ic348_simple_lc_scatter_brokenaxes,
    }
    cmap_dict = {"ngc": "jet", "ic": "jet_r"}

    for name in rs.names:

        followup_table = astropy.table.Table.read(
            os.path.join(lc_dir, name, f"Q0_objects_for_followup_{name}_2024_Oct_13.h5")
        )

        print(len(followup_table), " lightcurves incoming")

        for i, row in enumerate(followup_table):

            # for i, sid in enumerate(ds[bd & v0 & ~v_cand & ~per_bd].index):
            sid = row["SOURCEID"]
            S_value = row["var_Stetson_JHK"]
            Q_string = q_string(sid, rs.spread_dict[name], rs.q_dict[name])
            # print(f"sid={sid}, S={S_value:.2f}, Q={Q_string}")

            fig_lc = plot_dict[name](dat_dict[name], sid, cmap=cmap_dict[name])
            fig_lc.suptitle(
                f"{name}: {sid}, Q={Q_string}, S={S_value:.2f} \n Chisq J={row['var_J_red_chisq']:.2f} H={row['var_H_red_chisq']:.2f} K={row['var_K_red_chisq']:.2f}  \n Object name: {row['Name']}"
            )

            fig_lc.savefig(
                os.path.join(lc_dir, name, f"{i:03d}_{sid}_lc.png"),
                bbox_inches="tight",
            )


def generate_nonvariable_lightcurves_just_in_case():
    dat_dict = {
        "ngc": photometry.group_wserv_v2(
            photometry.load_wserv_v2(7, suffix="_outlier_cleaned_152")
        ),
        "ic": photometry.group_wserv_v2(photometry.load_wserv_v2(8)),
    }
    plot_dict = {
        "ngc": ngc1333_simple_lc_scatter_brokenaxes,
        "ic": ic348_simple_lc_scatter_brokenaxes,
    }
    cmap_dict = {"ngc": "jet", "ic": "jet_r"}

    for name in rs.names:
        match = rs.match_dict[name]

        nonvar_table = astropy.table.Table.read(
            os.path.join(
                non_lc_dir,
                name,
                f"nonvariable_objects_for_followup_{name}_2024_Oct_13.h5",
            )
        )

        # not_followup_table = match.approved[
        #     ~np.in1d(match.approved["SOURCEID"], followup_table["SOURCEID"])
        # ]

        print(len(nonvar_table), " lightcurves incoming")

        for i, row in enumerate(nonvar_table):

            # for i, sid in enumerate(ds[bd & v0 & ~v_cand & ~per_bd].index):
            sid = row["SOURCEID"]
            S_value = row["var_Stetson_JHK"]
            Q_string = q_string(sid, rs.spread_dict[name], rs.q_dict[name])
            # print(f"sid={sid}, S={S_value:.2f}, Q={Q_string}")

            fig_lc = plot_dict[name](dat_dict[name], sid, cmap=cmap_dict[name])
            fig_lc.suptitle(
                f"nonvar {name}: {sid}, Q={Q_string}, S={S_value:.2f} \n Chisq J={row['var_J_red_chisq']:.2f} H={row['var_H_red_chisq']:.2f} K={row['var_K_red_chisq']:.2f}  \n Object name: {row['Name']}"
            )

            fig_lc.savefig(
                os.path.join(non_lc_dir, name, f"{i:03d}_{sid}_lc.png"),
                bbox_inches="tight",
            )



if __name__ == "__main__":

    generate_list_for_followup()
    # generate_lightcurves_for_inspection()
    generate_nonvariable_lightcurves_just_in_case()
