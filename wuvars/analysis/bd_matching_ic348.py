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

### Part 1: loading the input catalogs.

aux_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/data/auxiliary_catalogs/IC348"

# Luhman et al. 2016, Table 1: Members of IC 348.
L16_T1_filepath = os.path.join(aux_path, "Luhman2016_Table1_apjaa2cabt1_mrt.txt")
L16_T1 = astropy.table.Table.read(L16_T1_filepath, format="ascii.cds")

# get the coordinates into a usable state
L16_T1_RA = Angle(
    (L16_T1["RAh"].data, L16_T1["RAm"].data, L16_T1["RAs"].data), unit=u.hourangle
)
L16_T1_Dec_sign = np.where(L16_T1["DE-"] == "+", 1, -1)
L16_T1_Dec = Angle(
    (L16_T1_Dec_sign * L16_T1["DEd"], L16_T1["DEm"], L16_T1["DEs"],), unit=u.deg,
)

L16_T1_coordinates = SkyCoord(ra=L16_T1_RA, dec=L16_T1_Dec)

# Part 2: loading our data.
# IC 348 is WSERV8
spread = spreadsheet.load_wserv_v2(8)
sm = spread["median"]
spreadsheet_coordinates = SkyCoord(
    ra=sm["RA"].values * u.rad, dec=sm["DEC"].values * u.rad
)

idx, d2d, d3d = L16_T1_coordinates.match_to_catalog_sky(spreadsheet_coordinates)


# Part 3: Doing the matching.

# Part 4: Outputting the results to file.

# SkyCoord(ra=(table['RAh'], table['RAm'], table['RAs']), dec=(table['DEd'], table['DEm'], table['DEs'])
