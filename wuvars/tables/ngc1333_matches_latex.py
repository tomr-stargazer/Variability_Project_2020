"""
Here is a script that takes in the brown dwarf matches in IC348 and makes a LaTeX table
from them.

"""

import os
import astropy.table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u

results_path = (
    "/Users/tsrice/Documents/Variability_Project_2020/Results/matched_catalogs/ngc1333"
)


bd_results_table = astropy.table.Table.read(
    os.path.join(results_path, "NGC1333_matches.dat"), format="ascii.ipac"
)

coords = SkyCoord(
    ra=bd_results_table["RA"] * u.rad, dec=bd_results_table["DEC"] * u.rad
)

# borrowed from https://stackoverflow.com/questions/16891340/remove-a-prefix-from-a-string
def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix) :]
    return text


names = bd_results_table["Name"]
for i in range(len(names)):
    names[i] = remove_prefix(names[i], "NGC 1333 ")


# we want our results table to have the following...
# NAME | SOURCEID | RA | Dec | J | H | K | SpT |

# for i, row in
