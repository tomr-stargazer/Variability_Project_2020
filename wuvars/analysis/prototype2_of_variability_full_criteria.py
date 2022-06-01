import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from wuvars.data import spreadsheet, photometry, quality_classes
import wuvars.analysis.variability_selection as sv
from wuvars.analysis.bd_matching_v2 import match_onc, match_ngc, match_ic
from wuvars.analysis.variability_selection_curved import (
    curve_Stetson,
    sv_jh,
    sv_hk,
    sv_jk,
    sv_jhk,
)

# general setup

spread = spreadsheet.load_v2()
stetson_directory = (
    "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/stetson_v_mag"
)
period_directory = (
    "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes"
)
subj_directory = period_directory
subj_sheets = pd.read_excel(
    os.path.join(subj_directory, "Q0_nonperiodic_followup.xlsx"), sheet_name=None
)

# in practice we don't use wserv 1 or 11, but I wrote this out and so I'll keep it
wserv_ids = [1, 5, 7, 8, 11]
n_obs_list = [130, 120, 171, 85, 110]

n_min_list = [60, 35, 80, 55, 65]
n_max_list = [100, 90, 160, 80, 100]

n_min_dict = {wserv: n_min for (wserv, n_min) in zip(wserv_ids, n_min_list)}
n_max_dict = {wserv: n_max for (wserv, n_max) in zip(wserv_ids, n_max_list)}

onc_match = match_onc()
ngc_match = match_ngc()
ic_match = match_ic()

reference_dict = {5: onc_match, 7: ngc_match, 8: ic_match}
short_names = {5: "ONC", 7: "NGC", 8: "IC"}
full_names = {5: "ONC", 7: "NGC 1333", 8: "IC 348"}


def select_targets(wserv, attr="approved"):

    sourceids = getattr(reference_dict[wserv], attr)["SOURCEID"]

    ds = spread[f"wserv{wserv}"]
    ref = np.in1d(ds.index, sourceids)

    return ref


def select_disks(wserv, attr="approved", choice="yes"):

    sourceids = getattr(reference_dict[wserv], attr)["SOURCEID"]
    disk_yn = getattr(reference_dict[wserv], attr)["IRexc"]

    disked_sourceids = sourceids[disk_yn == choice]

    ds = spread[f"wserv{wserv}"]
    ref = np.in1d(ds.index, disked_sourceids)

    return ref


def select_stetson_variables(wserv):

    n_min = n_min_dict[wserv]
    n_max = n_max_dict[wserv]

    stetson_result_997 = np.load(
        os.path.join(stetson_directory, f"WSERV{wserv}_result_grid_997.npy")
    )

    stetson_result_95 = np.load(
        os.path.join(stetson_directory, f"WSERV{wserv}_result_grid_95.npy")
    )

    ds = spread[f"wserv{wserv}"]
    q0 = sv.sq0(ds, n_min, n_max)
    #         q1 = sq1(ds, n_min, n_max)

    q1_j = sv.sq1_j(ds, n_min, n_max)
    q1_h = sv.sq1_h(ds, n_min, n_max)
    q1_k = sv.sq1_k(ds, n_min, n_max)

    q2 = sv.sq2(ds, n_min, n_max)

    suffix_997 = "_m_997"

    curve_Stetson(ds, stetson_result_997[0], stetson_result_997[1], "_m_997")
    curve_Stetson(ds, stetson_result_95[0], stetson_result_95[1], "_m_95")
    for bands in ["JH", "JK", "HK"]:
        curve_Stetson(
            ds,
            stetson_result_997[0],
            stetson_result_997[1],
            "_m_997",
            input_column=f"Stetson_{bands}",
        )

    v_jhk = sv_jhk(ds, Stetson_cutoff=0, suffix=suffix_997)
    #         v2 = sq2_variables()
    v_jh = sv_jh(ds, Stetson_cutoff=0, suffix=suffix_997)
    v_jk = sv_jk(ds, Stetson_cutoff=0, suffix=suffix_997)
    v_hk = sv_hk(ds, Stetson_cutoff=0, suffix=suffix_997)

    v2 = q2 & (v_jhk | v_jk | v_hk | v_jh)

    v1_jh = q1_j & q1_h & v_jh
    v1_jk = q1_j & q1_k & v_jk
    v1_hk = q1_h & q1_k & v_hk

    v1 = (
        (v1_jh & q1_j & q1_h)
        | (v1_hk & q1_k & q1_h)
        | (v1_jk & q1_j & q1_k)
        | (q2 & (v_jhk | v_jk | v_hk | v_jh))
    )

    return v2, v1


def select_periodic_variables(wserv):

    ds = spread[f"wserv{wserv}"]

    flags = ["Y", "Yw", "N", "YfY", "?fY", "YfYw", "?fYw", "YfN", "?fN"]
    periodic_flags = [flag for flag in flags if flag[-1] in ("Y", "w")]

    period_sheet = pd.read_excel(
        os.path.join(
            period_directory,
            f"{short_names[wserv]}_source_properties_periods_inspected.xlsx",
        )
    )
    periodics = period_sheet[np.in1d(period_sheet["Periodic?"], periodic_flags)]

    v_per = np.in1d(ds.index, periodics.values)

    return v_per


def select_subjective_variables(wserv):

    ds = spread[f"wserv{wserv}"]

    q0_variables = subj_sheets[short_names[wserv]]["SOURCEID"][
        subj_sheets[short_names[wserv]]["VARIABLE"] == 1
    ]
    # ngc_v0 = np.in1d(ngc_match.approved['SOURCEID'], q0_variables)
    v_subj = np.in1d(ds.index, q0_variables)

    return v_subj


def select_variables(wserv):

    ref = select_targets(wserv)

    v2, v1 = select_stetson_variables(wserv)
    v_per = select_periodic_variables(wserv)
    v_subjective = select_subjective_variables(wserv)

    vv = ref & (v2 | v1 | v_per | v_subjective)

    n_variables = np.sum(ref & (v2 | v1 | v_per | v_subjective))
    n_total = np.sum(ref)

    print(n_variables, n_total)

    return vv


if __name__ == "__main__":

    for wserv in [5, 7, 8]:
        print(f"WSERV{wserv}:")

        ref = select_targets(wserv)
        stat = select_targets(wserv, attr="statistical")
        v2, v1 = select_stetson_variables(wserv)
        v_per = select_periodic_variables(wserv)
        v_subj = select_subjective_variables(wserv)

        print(f"Ref: {np.sum(ref)} (Stat: {np.sum(stat)})")
        print(f"v2: {np.sum(ref & v2)} ({np.sum(stat & v2)})")
        print(f"v1: {np.sum(ref & v1)} ({np.sum(stat & v1)})")
        print(f"v_per: {np.sum(ref & v_per)} ({np.sum(stat & v_per)})")
        print(f"v_subj: {np.sum(ref & v_subj)} ({np.sum(stat & v_subj)})")

        vv = select_variables(wserv)

        print(f"v_total: {np.sum(ref & vv)} ({np.sum(stat & vv)})")
        rate = np.sum(stat & vv) / np.sum(stat)
        per_rate = np.sum(stat & v_per) / np.sum(stat)
        print(
            f"Statistical variability rate: {np.sum(stat & vv)}/{np.sum(stat)} = {rate:.2f}"
        )
        print(
            f"Statistical periodicity rate: {np.sum(stat & v_per)}/{np.sum(stat)} = {per_rate:.2f}"
        )

        print("")
        if wserv != 5:

            disks = select_disks(wserv)
            nondisks = select_disks(wserv, choice='no')    

            rate_disks = np.sum(stat & vv & disks) / np.sum(stat & disks)
            per_rate_disks = np.sum(stat & v_per & disks) / np.sum(stat & disks)
            rate_nondisks = np.sum(stat & vv & nondisks) / np.sum(stat & nondisks)
            per_rate_nondisks = np.sum(stat & v_per & nondisks) / np.sum(stat & nondisks)


            print(f"Disks: {np.sum(ref & disks)} (Stat: {np.sum(stat & disks)})")
            print(f"Non-disks: {np.sum(ref & nondisks)} (Stat: {np.sum(stat & nondisks)})")

            print(f"Variable Disks: {np.sum(ref & disks & vv)} (Stat: {np.sum(stat & disks & vv)})")
            print(f"Variable Non-disks: {np.sum(ref & nondisks & vv)} (Stat: {np.sum(stat & nondisks & vv)})")

            print(
                f"Statistical DISKED variability rate: {np.sum(stat & vv & disks)}/{np.sum(stat & disks)} = {rate_disks:.2f}"
            )
            print(
                f"Statistical DISKED periodicity rate: {np.sum(stat & v_per & disks)}/{np.sum(stat & disks)} = {per_rate_disks:.2f}"
            )
            print(
                f"Statistical NONDISKED variability rate: {np.sum(stat & vv & nondisks)}/{np.sum(stat & nondisks)} = {rate_nondisks:.2f}"
            )
            print(
                f"Statistical NONDISKED periodicity rate: {np.sum(stat & v_per & nondisks)}/{np.sum(stat & nondisks)} = {per_rate_nondisks:.2f}"
            )


        print("\n")
