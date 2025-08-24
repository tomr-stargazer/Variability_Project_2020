"""
Working on a unified bd_matching approach. This is v3.

Because the new paper only does NGC and IC, I'm only doing those here.

There are several different subsets which are defined in code below.
Each subset is a strict subset of the class immediately above it; 
imagine nesting Russian dolls, not a Venn diagram with complicated overlaps.

Here are my english-language definitions:
- all_matches: 
    every source which is
    (a) a confident member of the cluster in question*, and
    (b) matches within some matching radius** of a UKIRT source.
## WIP - still assessing what I'm calling these categories ##
## MLdwarfs ? MLmatches? Honestly, it doesn't matter; I'll just use approved
## approved_v3 ... or just approved
## statistical
- M4:
    subset of the above which has SpT >= M4.0 .
- approved:
    subset of the above which
    (a) upon my visual inspection of its lightcurve, did not have data 
        poor enough to say "there is no useful variability information here"
- statistical:
    subset of the above which
    (a) has Q1_k quality (at least) and
    (b) has median K error below 0.05 mag
- color:
    subset of the above which
    (a) has Q2 quality

Our goal:
- check for periods and "interesting" lightcurves among all `approved`
- do K-band variability amplitude stats on `statistical`
- search for color variability stats among `color` 

That's the paper. Periods, interesting lightcurves, K band amplitude stats, and 
color variability.

* 'confident member' means:
(a) listed in Tables 1 or 2 of Luhman+16 (for Perseus), or
(b) Bayes Factor > 99 i.e., ["log(BF)"] > np.log10(99) (for Orion)

** matching radius is 0.5'', chosen by histogram analysis (by eye).
This supersedes the 0.37'' radius I used in v2.

"""


import os

import astropy.table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from wuvars.analysis.luhman16_coord_handler import coords_from_Luhman_table
from wuvars.analysis.spectral_type_to_number import (get_num_from_SpT,
                                                     get_SpT_from_num)
from wuvars.analysis.spectral_type_to_temperature import (get_SpT_from_Teff,
                                                          get_Teff_from_SpT)
from wuvars.data import photometry, quality_classes, spreadsheet


class MatchStruct:
    pass

# The following is the "old" approved list
lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"
inspect_ngc = pd.read_csv(
    os.path.join(lc_dir, "inspection_ngc.csv"), skipinitialspace=True
)
inspect_ic = pd.read_csv(
    os.path.join(lc_dir, "inspection_ic.csv"), skipinitialspace=True
)

rejected_sources_ngc_old = inspect_ngc["SOURCEID"][inspect_ngc["exclude?"] == "yes"]
rejected_sources_ic_old = inspect_ic["SOURCEID"][inspect_ic["exclude?"] == "yes"]

approved_sources_ngc_old = inspect_ngc["SOURCEID"][inspect_ngc["exclude?"] != "yes"]
approved_sources_ic_old = inspect_ic["SOURCEID"][inspect_ic["exclude?"] != "yes"]

# The following is the "new" approved list
lc_dir_v3 = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3.5"
inspect_ngc_v3 = pd.read_csv(
    os.path.join(lc_dir_v3, "new_inspection_ngc.csv"), skipinitialspace=True
)
inspect_ic_v3 = pd.read_csv(
    os.path.join(lc_dir_v3, "new_inspection_ic.csv"), skipinitialspace=True
)

rejected_sources_ngc_new = inspect_ngc_v3["SOURCEID"][inspect_ngc_v3["exclude?"] == "yes"]
rejected_sources_ic_new = inspect_ic_v3["SOURCEID"][inspect_ic_v3["exclude?"] == "yes"]

approved_sources_ngc_new = inspect_ngc_v3["SOURCEID"][inspect_ngc_v3["exclude?"] != "yes"]
approved_sources_ic_new = inspect_ic_v3["SOURCEID"][inspect_ic_v3["exclude?"] != "yes"]

approved_sources_ngc = pd.concat([approved_sources_ngc_new, approved_sources_ngc_old])
approved_sources_ic = pd.concat([approved_sources_ic_new, approved_sources_ic_old])

def match_ngc():

    ngc_match = MatchStruct()

    # =========================================== #
    # === Part 1: loading the input catalogs. === #
    # =========================================== #

    aux_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/data/auxiliary_catalogs/NGC1333"

    # Luhman et al. 2016, Table 2: Members of NGC 1333.
    L16_T2_filepath = os.path.join(aux_path, "Luhman2016_Table2_apjaa2cabt1_mrt.txt")
    L16_T2 = astropy.table.Table.read(L16_T2_filepath, format="ascii.cds")

    # get the coordinates into a usable state
    L16_T2_coordinates = coords_from_Luhman_table(L16_T2)

    L16_T2_index = np.arange(len(L16_T2))
    L16_T2.add_column(L16_T2_index, index=0, name="L16_T2_index")

    # This column contains the 'adopted' spectral type for each source.
    L16_SpT = L16_T2["Adopt"]

    # onc_spt = np.array([get_SpT_from_Teff(x) for x in lowmass_matched["Teff"]])
    # not_lowmass_onc_spt = np.array([get_SpT_from_Teff(x) for x in matched['Teff']])
    L16_SpT_num = []
    for SpT in L16_SpT:
        try:
            L16_SpT_num.append(get_num_from_SpT(SpT))
        # This exception applies to any 'masked' entries, which should be treated as NaN.
        except IndexError:
            L16_SpT_num.append(np.nan)


    L16_T2["SpT"] = L16_SpT_num
    ngc_Teff = np.array([get_Teff_from_SpT(x) for x in L16_T2["SpT"]])

    L16_T2["Teff"] = ngc_Teff
    L16_T2["RA"] = np.radians(L16_T2_coordinates.ra)
    L16_T2["DEC"] = np.radians(L16_T2_coordinates.dec)

    ngc_match.input_catalog = L16_T2

    # ======================================= #
    # === Part 2: loading our UKIRT data. === #
    # ======================================= #

    # NGC 1333 is WSERV7
    spread = spreadsheet.load_wserv_v2(7)
    sm = spread["median"]
    sv = spread["variability"]
    ss = spread["std"]
    sr = spread["range_9010"]
    spreadsheet_coordinates = SkyCoord(
        ra=sm["RA"].values * u.rad, dec=sm["DEC"].values * u.rad
    )

    ngc_q = quality_classes.load_q(7)

    # =================================== #
    # === Part 3: Doing the matching. === #
    # =================================== #

    idx, d2d, d3d = L16_T2_coordinates.match_to_catalog_sky(spreadsheet_coordinates)

    # Maximum separation of 0.5'' chosen by inspecting a histogram of closest matches.
    max_sep = 0.5 * u.arcsec
    sep_constraint = d2d < max_sep

    matches_sm = sm.iloc[idx[sep_constraint]]
    matches_sm = matches_sm.rename(columns=lambda name: "median_" + name)
    matches_sv = sv.iloc[idx[sep_constraint]]
    matches_sv = matches_sv.rename(columns=lambda name: "var_" + name)
    matches_ss = ss.iloc[idx[sep_constraint]]
    matches_ss = matches_ss.rename(columns=lambda name: "std_" + name)
    matches_sr = sr.iloc[idx[sep_constraint]]
    matches_sr = matches_sr.rename(columns=lambda name: "range_" + name)

    matched = L16_T2[sep_constraint]
    unmatched = L16_T2[~sep_constraint]

    joint_matches = astropy.table.hstack(
        [
            astropy.table.Table.from_pandas(matches_sm),
            astropy.table.Table.from_pandas(matches_sv),
            astropy.table.Table.from_pandas(matches_ss),
            astropy.table.Table.from_pandas(matches_sr),
            matched,
        ]
    )
    joint_matches.add_column(matches_sm.index, index=0, name="SOURCEID")
    joint_matches.add_column(d2d[sep_constraint].to(u.arcsec), name="d2d_arcsec")

    ngc_match.all_matches = joint_matches

    lowmass_criteria = joint_matches["SpT"] >= 4.5
    ngc_match.lowmass = joint_matches[lowmass_criteria]
    ngc_match.not_lowmass = joint_matches[~lowmass_criteria]

    # New as of Oct 2023.
    M4_criteria = joint_matches["SpT"] >= 4.0
    ngc_match.M4 = joint_matches[M4_criteria]
    ngc_match.not_M4 = joint_matches[~M4_criteria]

    ML_criteria = joint_matches["SpT"] >= 0.0
    ngc_match.ML = joint_matches[ML_criteria]
    ngc_match.not_ML = joint_matches[~ML_criteria]

    approved_indices_ngc = np.in1d(joint_matches["SOURCEID"], approved_sources_ngc)

    approved_criteria = ML_criteria & approved_indices_ngc
    rejected_criteria = ML_criteria & ~approved_indices_ngc

    ngc_match.approved = joint_matches[approved_criteria]
    ngc_match.rejected = joint_matches[rejected_criteria]

    statistical_criteria = (
        approved_criteria
        & (joint_matches["median_KAPERMAG3ERR"] < 0.05)
        & ngc_q.q1_k[joint_matches["SOURCEID"]].values
    )
    ngc_match.statistical = joint_matches[statistical_criteria]

    color_criteria = statistical_criteria & ngc_q.q2[joint_matches["SOURCEID"]].values
    ngc_match.color = joint_matches[color_criteria]

    ngc_match.unmatched = L16_T2[~sep_constraint & (L16_T2["SpT"] >= 0.0)]

    return ngc_match


def match_ic():

    ic_match = MatchStruct()

    # =========================================== #
    # === Part 1: loading the input catalogs. === #
    # =========================================== #

    aux_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/data/auxiliary_catalogs/IC348"

    # Luhman et al. 2016, Table 1: Members of IC 348.
    L16_T1_filepath = os.path.join(aux_path, "Luhman2016_Table1_apjaa2cabt1_mrt.txt")
    L16_T1 = astropy.table.Table.read(L16_T1_filepath, format="ascii.cds")

    # get the coordinates into a usable state
    L16_T1_coordinates = coords_from_Luhman_table(L16_T1)

    L16_T1_index = np.arange(len(L16_T1))
    L16_T1.add_column(L16_T1_index, index=0, name="L16_T1_index")

    # This column contains the 'adopted' spectral type for each source.
    L16_SpT = L16_T1["Adopt"]
    L16_SpT_num = []
    for SpT in L16_SpT:
        try:
            L16_SpT_num.append(get_num_from_SpT(SpT))
        # This exception applies to any 'masked' entries, which should be treated as NaN.
        except IndexError:
            L16_SpT_num.append(np.nan)

    L16_T1["SpT"] = L16_SpT_num
    ic_Teff = np.array([get_Teff_from_SpT(x) for x in L16_T1["SpT"]])

    L16_T1["Teff"] = ic_Teff
    L16_T1["RA"] = np.radians(L16_T1_coordinates.ra)
    L16_T1["DEC"] = np.radians(L16_T1_coordinates.dec)

    ic_match.input_catalog = L16_T1
    
    # ======================================= #
    # === Part 2: loading our UKIRT data. === #
    # ======================================= #

    # IC 348 is WSERV8
    spread = spreadsheet.load_wserv_v2(8)
    sm = spread["median"]
    sv = spread["variability"]
    ss = spread["std"]
    sr = spread["range_9010"]
    spreadsheet_coordinates = SkyCoord(
        ra=sm["RA"].values * u.rad, dec=sm["DEC"].values * u.rad
    )

    ic_q = quality_classes.load_q(8)

    # =================================== #
    # === Part 3: Doing the matching. === #
    # =================================== #

    idx, d2d, d3d = L16_T1_coordinates.match_to_catalog_sky(spreadsheet_coordinates)

    # Maximum separation of 0.5'' chosen by inspecting a histogram of closest matches.
    max_sep = 0.5 * u.arcsec
    sep_constraint = d2d < max_sep

    matches_sm = sm.iloc[idx[sep_constraint]]
    matches_sm = matches_sm.rename(columns=lambda name: "median_" + name)
    matches_sv = sv.iloc[idx[sep_constraint]]
    matches_sv = matches_sv.rename(columns=lambda name: "var_" + name)
    matches_ss = ss.iloc[idx[sep_constraint]]
    matches_ss = matches_ss.rename(columns=lambda name: "std_" + name)
    matches_sr = sr.iloc[idx[sep_constraint]]
    matches_sr = matches_sr.rename(columns=lambda name: "range_" + name)

    matched = L16_T1[sep_constraint]
    unmatched = L16_T1[~sep_constraint]
    joint_matches = astropy.table.hstack(
        [
            astropy.table.Table.from_pandas(matches_sm),
            astropy.table.Table.from_pandas(matches_sv),
            astropy.table.Table.from_pandas(matches_ss),
            astropy.table.Table.from_pandas(matches_sr),
            matched,
        ]
    )
    joint_matches.add_column(matches_sm.index, index=0, name="SOURCEID")
    joint_matches.add_column(d2d[sep_constraint].to(u.arcsec), name="d2d_arcsec")

    ic_match.all_matches = joint_matches

    lowmass_criteria = joint_matches["SpT"] >= 4.5
    ic_match.lowmass = joint_matches[lowmass_criteria]
    ic_match.not_lowmass = joint_matches[~lowmass_criteria]

    # New as of Oct 2023.
    M4_criteria = joint_matches["SpT"] >= 4.0
    ic_match.M4 = joint_matches[M4_criteria]
    ic_match.not_M4 = joint_matches[~M4_criteria]

    ML_criteria = joint_matches["SpT"] >= 0.0
    ic_match.ML = joint_matches[ML_criteria]
    ic_match.not_ML = joint_matches[~ML_criteria]

    approved_indices_ic = np.in1d(joint_matches["SOURCEID"], approved_sources_ic)

    approved_criteria = ML_criteria & approved_indices_ic
    rejected_criteria = ML_criteria & ~approved_indices_ic

    ic_match.approved = joint_matches[approved_criteria]
    ic_match.rejected = joint_matches[rejected_criteria]

    statistical_criteria = (
        approved_criteria
        & (joint_matches["median_KAPERMAG3ERR"] < 0.05)
        & ic_q.q1_k[joint_matches["SOURCEID"]].values
    )
    ic_match.statistical = joint_matches[statistical_criteria]

    color_criteria = statistical_criteria & ic_q.q2[joint_matches["SOURCEID"]].values
    ic_match.color = joint_matches[color_criteria]

    ic_match.unmatched = L16_T1[~sep_constraint & (L16_T1["SpT"] >= 0.0)]

    return ic_match
