"""
This is a script which cross-matches all relevant input BD catalogs 
in IC 348 to our data.

I'm leaning on the prototype notebook "Playing with Luhman et al. 2016"
for much of this.

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

# Choose all objects with spectral type later than M6.0
L16_bds = np.array(L16_SpT_num) >= 6.0
L16_lowmass_ic = (np.array(L16_SpT_num) >= 4.5) 


# ======================================= #
# === Part 2: loading our UKIRT data. === #
# ======================================= #

# IC 348 is WSERV8
spread = spreadsheet.load_wserv_v2(8)
sm = spread["median"]
spreadsheet_coordinates = SkyCoord(
    ra=sm["RA"].values * u.rad, dec=sm["DEC"].values * u.rad,
)


# =================================== #
# === Part 3: Doing the matching. === #
# =================================== #

idx, d2d, d3d = L16_T1_coordinates.match_to_catalog_sky(spreadsheet_coordinates)

# Maximum separation of 0.37'' chosen by inspecting a histogram of closest matches.
max_sep = 0.37 * u.arcsec
sep_constraint = d2d < max_sep

# We're going to compute
# (a) all matched IC348 members, for quality control / comparisons, and
# (b) just the brown dwarfs

matches = sm.iloc[idx[sep_constraint]]
matched = L16_T1[sep_constraint]
joint_matches = astropy.table.hstack(
    [astropy.table.Table.from_pandas(matches), matched]
)
joint_matches.add_column(matches.index, index=0, name="SOURCEID")

bd_matches = sm.iloc[idx[sep_constraint & L16_bds]]
bd_matched = L16_T1[sep_constraint & L16_bds]
bd_joint_matches = astropy.table.hstack(
    [astropy.table.Table.from_pandas(bd_matches), bd_matched]
)
bd_joint_matches.add_column(bd_matches.index, index=0, name="SOURCEID")

lowmass_ic_matches = sm.iloc[idx[sep_constraint & L16_lowmass_ic]]
lowmass_ic_matched = L16_T1[sep_constraint & L16_lowmass_ic]
lowmass_ic_joint_matches = astropy.table.hstack(
    [astropy.table.Table.from_pandas(lowmass_ic_matches), lowmass_ic_matched]
)
lowmass_ic_joint_matches.add_column(lowmass_ic_matches.index, index=0, name="SOURCEID")


# =============================================== #
# === Part 4: Outputting the results to file. === #
# =============================================== #

results_path = (
    "/Users/tsrice/Documents/Variability_Project_2020/Results/matched_catalogs/ic348"
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
    "Name", "SOURCEID", "RA", "DEC", j, h, k, jerr, herr, kerr, "Adopt"
]
bd_results_table["Adopt"].name = "SpT"

# we want our results table to have the following...
# NAME | SOURCEID | RA | Dec | J | H | K | SpT |

bd_results_table.write(
    os.path.join(results_path, "IC348_matches.dat"), format="ascii.ipac", overwrite=True
)

# SkyCoord(ra=(table['RAh'], table['RAm'], table['RAs']), dec=(table['DEd'], table['DEm'], table['DEs'])
