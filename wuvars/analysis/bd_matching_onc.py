"""
This is a script which cross-matches all relevant input BD catalogs 
in the ONC to our data.

I'm leaning on the script `bd_matching_ic348.py` for much of this.

At the moment I am prototyping this on using only Robberto et al. 2020's catalog, 
which might not be totally optimal.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u

from wuvars.data import spreadsheet, photometry
from wuvars.analysis.spectral_type_to_number import get_num_from_SpT
from wuvars.analysis.luhman16_coord_handler import coords_from_Luhman_table


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


# # Luhman et al. 2016, Table 2: Members of NGC 1333.
# L16_T2_filepath = os.path.join(aux_path, "Luhman2016_Table2_apjaa2cabt1_mrt.txt")
# L16_T2 = astropy.table.Table.read(L16_T2_filepath, format="ascii.cds")

# # get the coordinates into a usable state
# L16_T2_coordinates = coords_from_Luhman_table(L16_T2)
robb_coordinates = SkyCoord(ra=robb_joined["RAdeg"], dec=robb_joined["DEdeg"])

# # This column contains the 'adopted' spectral type for each source.
# L16_SpT = L16_T2["Adopt"]
# L16_SpT_num = []
# for SpT in L16_SpT:
#     try:
#         L16_SpT_num.append(get_num_from_SpT(SpT))
#     # This exception applies to any 'masked' entries, which should be treated as NaN.
#     except IndexError:
#         L16_SpT_num.append(np.nan)

# # Choose all objects with spectral type later than M6.0
# L16_bds = np.array(L16_SpT_num) >= 6.0
confident_members = robb_joined["log(BF)"] > np.log10(99)
confident_bds = (robb_joined["Mstar"] < 0.08) & (robb_joined["log(BF)"] > np.log10(99))

confident_lowmass = (robb_joined["Teff"] <= 3200) & (robb_joined["log(BF)"] > np.log10(99))


# ======================================= #
# === Part 2: loading our UKIRT data. === #
# ======================================= #

# ONC is WSERV5
spread = spreadsheet.load_wserv_v2(5)
sm = spread["median"]
spreadsheet_coordinates = SkyCoord(
    ra=sm["RA"].values * u.rad, dec=sm["DEC"].values * u.rad
)


# =================================== #
# === Part 3: Doing the matching. === #
# =================================== #

idx, d2d, d3d = robb_coordinates.match_to_catalog_sky(spreadsheet_coordinates)


# Maximum separation of 0.37'' chosen by inspecting a histogram of closest matches.
max_sep = 0.37 * u.arcsec
sep_constraint = d2d < max_sep

# # We're going to compute
# # (a) all matched NGC1333 members, for quality control / comparisons, and
# # (b) just the brown dwarfs

matches = sm.iloc[idx[sep_constraint & confident_members]]
matched = robb_joined[sep_constraint & confident_members]
joint_matches = astropy.table.hstack(
    [astropy.table.Table.from_pandas(matches), matched]
)
joint_matches.add_column(matches.index, index=0, name="SOURCEID")

bd_matches = sm.iloc[idx[sep_constraint & confident_bds]]
bd_matched = robb_joined[sep_constraint & confident_bds]
bd_joint_matches = astropy.table.hstack(
    [astropy.table.Table.from_pandas(bd_matches), bd_matched]
)
bd_joint_matches.add_column(bd_matches.index, index=0, name="SOURCEID")


lowmass_matches = sm.iloc[idx[sep_constraint & confident_lowmass]]
lowmass_matches_spread = spread.iloc[idx[sep_constraint & confident_lowmass]]
lowmass_matched = robb_joined[sep_constraint & confident_lowmass]
lowmass_joint_matches = astropy.table.hstack(
    [astropy.table.Table.from_pandas(lowmass_matches), lowmass_matched]
)
lowmass_joint_matches.add_column(lowmass_matches.index, index=0, name="SOURCEID")

# # =============================================== #
# # === Part 4: Outputting the results to file. === #
# # =============================================== #

results_path = (
    "/Users/tsrice/Documents/Variability_Project_2020/Results/matched_catalogs/onc"
)

j = "JAPERMAG3"
h = "HAPERMAG3"
k = "KAPERMAG3"
jmh = "JMHPNT"
hmk = "HMKPNT"

jerr = j + "ERR"
herr = h + "ERR"
kerr = k + "ERR"

bd_results_table = bd_joint_matches[
    "index", "SOURCEID", "RA", "DEC", j, h, k, jerr, herr, kerr, "Mstar", "Teff"
]
# bd_results_table["Adopt"].name = "SpT"

# we want our results table to have the following...
# NAME | SOURCEID | RA | Dec | J | H | K | Mass | Teff

bd_results_table.write(
    os.path.join(results_path, "ONC_matches.dat"),
    format="ascii.ipac",
    overwrite=True,
)

# # SkyCoord(ra=(table['RAh'], table['RAm'], table['RAs']), dec=(table['DEd'], table['DEm'], table['DEs'])
