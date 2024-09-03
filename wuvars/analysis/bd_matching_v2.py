"""
Working on a unified bd_matching approach.

Prototyping it in script form here.

There are several different subsets which are defined in code below.
Each subset is a strict subset of the class immediately above it; 
imagine nesting Russian dolls, not a Venn diagram with complicated overlaps.

Here are my english-language definitions:
- all_matches: 
    every source which is
    (a) a confident member of the cluster in question*, and
    (b) matches within some matching radius** of a UKIRT source.
- lowmass:
    subset of the above which has either
    (a) Teff < 3200 (if Teff is primary), or 
    (b) SpT >= M4.5 (if SpT is primary).
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

** matching radius is 0.37'', chosen by histogram analysis (by eye).

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


lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"
inspect_onc = pd.read_csv(
    os.path.join(lc_dir, "inspection_onc.csv"), skipinitialspace=True
)
inspect_ngc = pd.read_csv(
    os.path.join(lc_dir, "inspection_ngc.csv"), skipinitialspace=True
)
inspect_ic = pd.read_csv(
    os.path.join(lc_dir, "inspection_ic.csv"), skipinitialspace=True
)

rejected_sources_onc = inspect_onc["SOURCEID"][inspect_onc["exclude?"] == "yes"]
rejected_sources_ngc = inspect_ngc["SOURCEID"][inspect_ngc["exclude?"] == "yes"]
rejected_sources_ic = inspect_ic["SOURCEID"][inspect_ic["exclude?"] == "yes"]

approved_sources_onc = inspect_onc["SOURCEID"][inspect_onc["exclude?"] != "yes"]
approved_sources_ngc = inspect_ngc["SOURCEID"][inspect_ngc["exclude?"] != "yes"]
approved_sources_ic = inspect_ic["SOURCEID"][inspect_ic["exclude?"] != "yes"]


def match_onc():

    onc_match = MatchStruct()

    # =========================================== #
    # === Part 1: loading the input catalogs. === #
    # =========================================== #

    aux_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/data/auxiliary_catalogs/ONC"

    # Robberto et al. 2020, Table 2: WFC3-IR Photometry of ONC Sources
    filepath_table2 = os.path.join(aux_path, "Robberto_2020_Table2_apjab911et2_mrt.txt")
    robb_table2 = astropy.table.Table.read(filepath_table2, format="ascii")[:-1]

    # Robberto et al. 2020, Table 4: Estimated Parameters of the ONC Sample
    filepath_table4 = os.path.join(aux_path, "Robberto_2020_Table4_apjab911et4_mrt.txt")
    robb_table4 = astropy.table.Table.read(filepath_table4, format="ascii")

    robb_joined = astropy.table.join(robb_table2, robb_table4)
    onc_spt = np.array([get_SpT_from_Teff(x) for x in robb_joined["Teff"]])
    robb_joined["SpT"] = onc_spt

    robb_coordinates = SkyCoord(ra=robb_joined["RAdeg"], dec=robb_joined["DEdeg"])

    # # This column contains the Teff.
    confident_members = robb_joined["log(BF)"] > np.log10(99)

    # ======================================= #
    # === Part 2: loading our UKIRT data. === #
    # ======================================= #

    # ONC is WSERV5
    spread = spreadsheet.load_wserv_v2(5)
    sm = spread["median"]
    sv = spread["variability"]
    ss = spread["std"]
    sr = spread["range_9010"]
    spreadsheet_coordinates = SkyCoord(
        ra=sm["RA"].values * u.rad, dec=sm["DEC"].values * u.rad
    )

    onc_q = quality_classes.load_q(5)

    # =================================== #
    # === Part 3: Doing the matching. === #
    # =================================== #

    idx, d2d, d3d = robb_coordinates.match_to_catalog_sky(spreadsheet_coordinates)

    # Maximum separation of 0.37'' chosen by inspecting a histogram of closest matches.
    max_sep = 0.37 * u.arcsec
    sep_constraint = d2d < max_sep

    matches_sm = sm.iloc[idx[sep_constraint & confident_members]]
    matches_sm = matches_sm.rename(columns=lambda name: "median_" + name)
    matches_sv = sv.iloc[idx[sep_constraint & confident_members]]
    matches_sv = matches_sv.rename(columns=lambda name: "var_" + name)
    matches_ss = ss.iloc[idx[sep_constraint & confident_members]]
    matches_ss = matches_ss.rename(columns=lambda name: "std_" + name)
    matches_sr = sr.iloc[idx[sep_constraint & confident_members]]
    matches_sr = matches_sr.rename(columns=lambda name: "range_" + name)

    matched = robb_joined[sep_constraint & confident_members]
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

    onc_match.all_matches = joint_matches

    lowmass_criteria = joint_matches["Teff"] < 3200
    onc_match.lowmass = joint_matches[lowmass_criteria]
    onc_match.not_lowmass = joint_matches[~lowmass_criteria]

    approved_indices_onc = np.in1d(joint_matches["SOURCEID"], approved_sources_onc)

    approved_criteria = lowmass_criteria & approved_indices_onc
    onc_match.approved = joint_matches[approved_criteria]

    statistical_criteria = (
        approved_criteria
        & (joint_matches["median_KAPERMAG3ERR"] < 0.05)
        & onc_q.q1_k[joint_matches["SOURCEID"]].values
    )
    onc_match.statistical = joint_matches[statistical_criteria]

    color_criteria = statistical_criteria & onc_q.q2[joint_matches["SOURCEID"]].values
    onc_match.color = joint_matches[color_criteria]

    return onc_match


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

    # Maximum separation of 0.37'' chosen by inspecting a histogram of closest matches.
    max_sep = 0.37 * u.arcsec
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

    ngc_match.all_matches = joint_matches

    lowmass_criteria = joint_matches["SpT"] >= 4.5
    ngc_match.lowmass = joint_matches[lowmass_criteria]
    ngc_match.not_lowmass = joint_matches[~lowmass_criteria]

    # New as of Oct 2023.
    M4_criteria = joint_matches["SpT"] >= 4.0
    ngc_match.M4 = joint_matches[M4_criteria]
    ngc_match.not_M4 = joint_matches[~M4_criteria]

    approved_indices_ngc = np.in1d(joint_matches["SOURCEID"], approved_sources_ngc)

    approved_criteria = lowmass_criteria & approved_indices_ngc
    ngc_match.approved = joint_matches[approved_criteria]

    statistical_criteria = (
        approved_criteria
        & (joint_matches["median_KAPERMAG3ERR"] < 0.05)
        & ngc_q.q1_k[joint_matches["SOURCEID"]].values
    )
    ngc_match.statistical = joint_matches[statistical_criteria]

    color_criteria = statistical_criteria & ngc_q.q2[joint_matches["SOURCEID"]].values
    ngc_match.color = joint_matches[color_criteria]

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

    # Maximum separation of 0.37'' chosen by inspecting a histogram of closest matches.
    max_sep = 0.37 * u.arcsec
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

    ic_match.all_matches = joint_matches

    lowmass_criteria = joint_matches["SpT"] >= 4.5
    ic_match.lowmass = joint_matches[lowmass_criteria]
    ic_match.not_lowmass = joint_matches[~lowmass_criteria]

    # New as of Oct 2023.
    M4_criteria = joint_matches["SpT"] >= 4.0
    ic_match.M4 = joint_matches[M4_criteria]
    ic_match.not_M4 = joint_matches[~M4_criteria]

    approved_indices_ic = np.in1d(joint_matches["SOURCEID"], approved_sources_ic)

    approved_criteria = lowmass_criteria & approved_indices_ic
    ic_match.approved = joint_matches[approved_criteria]

    statistical_criteria = (
        approved_criteria
        & (joint_matches["median_KAPERMAG3ERR"] < 0.05)
        & ic_q.q1_k[joint_matches["SOURCEID"]].values
    )
    ic_match.statistical = joint_matches[statistical_criteria]

    color_criteria = statistical_criteria & ic_q.q2[joint_matches["SOURCEID"]].values
    ic_match.color = joint_matches[color_criteria]

    return ic_match
