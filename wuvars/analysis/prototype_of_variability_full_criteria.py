# variable selection - a script.

# imagine we are loading things up just for NGC 1333.

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wuvars.analysis.variability_selection as sv
from wuvars.analysis.bd_matching_v2 import match_ngc
from wuvars.analysis.variability_selection_curved import (curve_Stetson, sv_hk,
                                                          sv_jh, sv_jhk, sv_jk)
from wuvars.data import photometry, quality_classes, spreadsheet

spread = spreadsheet.load_v2()

ngc_match = match_ngc()

load_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/stetson_v_mag"

StetsonJHK_dict = {}
StetsonJHK_dict_95 = {}
StetsonJHK_dict_50 = {}
wserv = 7
r997 = np.load(
    os.path.join(load_dir, f"WSERV{wserv}_result_grid_997.npy")
)
StetsonJHK_dict[wserv] = r997
r95 = np.load(
    os.path.join(load_dir, f"WSERV{wserv}_result_grid_95.npy")
)
StetsonJHK_dict_95[wserv] = r95    
r50 = np.load(
    os.path.join(load_dir, f"WSERV{wserv}_result_grid_50.npy")
)
StetsonJHK_dict_50[wserv] = r50     

flags = ['Y', 'Yw', 'N', 'YfY', '?fY', 'YfYw', '?fYw', "YfN", "?fN"]
periodic_flags = [flag for flag in flags if flag[-1] in ('Y', 'w')]
spreadsheet_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes"

period_sheet_ngc = pd.read_excel(
    os.path.join(spreadsheet_dir, "NGC_source_properties_periods_inspected.xlsx")
)
periodic_ngc = period_sheet_ngc[np.in1d(period_sheet_ngc['Periodic?'], periodic_flags)]

periodic = {}
periodic[7] = periodic_ngc

lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"

inspect_ngc = pd.read_csv(os.path.join(lc_dir, "inspection_ngc.csv"), skipinitialspace=True)

approved_sources_ngc = inspect_ngc['SOURCEID'][inspect_ngc['exclude?'] != 'yes']
approved = {}
approved[7] = approved_sources_ngc

# these are some parameters that are considered "input"
wserv = 7
n_min = 80
n_max = 160


ds = getattr(spread, f"wserv{wserv}")

q0 = sv.sq0(ds, n_min, n_max)
#         q1 = sq1(ds, n_min, n_max)

q1_j = sv.sq1_j(ds, n_min, n_max)
q1_h = sv.sq1_h(ds, n_min, n_max)
q1_k = sv.sq1_k(ds, n_min, n_max)

q2 = sv.sq2(ds, n_min, n_max)

#         v0 = sq0_variables(ds, n_min, n_max, Stetson_cutoff=S)
#         v1 = sq1_variables(ds, n_min, n_max, Stetson_cutoff=S)
#         v2 = sq2_variables(ds, n_min, n_max, Stetson_cutoff=S)
suffix_997 = "_m_997"
suffix_95 = "_m_95"
curve_Stetson(ds, StetsonJHK_dict[wserv][0], StetsonJHK_dict[wserv][1], "_m_997")
curve_Stetson(ds, StetsonJHK_dict_95[wserv][0], StetsonJHK_dict_95[wserv][1], "_m_95")
for bands in ["JH", "JK", "HK"]:
    curve_Stetson(
        ds,
        StetsonJHK_dict[wserv][0],
        StetsonJHK_dict[wserv][1],
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

v_ = v1 | v2

v0 = (v_jhk | v1_jh | v1_jk | v1_hk) & ~v1 & ~v2

v_jhk_cand = sv_jhk(ds, Stetson_cutoff=0, suffix=suffix_95)

v_cand = q2 & v_jhk_cand & ~v_

bd = np.in1d(ds.index, approved[wserv].values)
per_bd = np.in1d(ds.index, periodic[wserv].values)

print(f"WSERV{wserv}:")
print(f"  VLMS stars in our data: {np.sum(bd)}")
print(f"  VLMS Q=2 stars: {np.sum(bd & q2)}")
print(f"  VLMS Q=2 variables (ignoring periodicity): {np.sum(bd & q2 & v2)}")
print(f"  VLMS Q=2 'almost' variables: {np.sum(bd & q2 & v_jhk_cand & ~v2)}")
print(f"  VLMS Q=2 periodic: {np.sum(bd & q2 & per_bd)}")
print(f"  VLMS Q=2 variables (incl. periodics): {np.sum(bd & q2 & (v2 | per_bd))}")
print("")
print(f"BDs:                {np.sum(bd)}")
print(f"BDs (Q=2):          {np.sum(bd & q2)}")
print(f"Periodic BDs:       {np.sum(bd & per_bd)}")
print(f"Periodic BDs (Q=2): {np.sum(bd & per_bd & q2)}")

print(f"Variable BDs (Q=2): {np.sum(bd & v2)}")
print(f"  '' - periodics  : {np.sum(bd & v2 & ~per_bd)}")
print(f"Variable BDs (Q=1): {np.sum(bd & v1)}")
print(f"  '' - v2:          {np.sum(bd & v1 & ~v2)}")

print(f"V0 BDs (Q=0):       {np.sum(bd & v0)}")
print(f"Variable BDs (Q=1+2+per): {np.sum(bd & (v1 | v2 | per_bd))}")
print("\n\n")


print("NGC 1333:")
for attr, value in ngc_match.__dict__.items():
    if attr != 'not_lowmass':
        print(f"{attr:13s}", len(value))

sheets = pd.read_excel("Q0_nonperiodic_followup.xlsx", sheet_name=None)

q0_variables = sheets['NGC']['SOURCEID'][sheets['NGC']['VARIABLE'] == 1]

ngc_v0 = np.in1d(ngc_match.approved['SOURCEID'], q0_variables)
print(len(ngc_match.approved[ngc_v0]))        
v_per = np.in1d(ds.index, periodic[wserv].values)

v0_conf = np.in1d(ds.index, q0_variables)

ref = np.in1d(ds.index, ngc_match.approved['SOURCEID'])
print(np.sum(ref & q2), np.sum(ref & v2), np.sum(ref & v1), np.sum(ref & v_per), np.sum(ref & q2 & v_per), np.sum(ref & v0_conf))


print(np.sum(ref & (v2 | v1 | v_per | v0_conf)), np.sum(ref))
