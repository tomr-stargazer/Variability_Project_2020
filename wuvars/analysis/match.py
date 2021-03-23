"""
Given a coordinate pair, checks in 

See also this astropy documentation guide: https://learn.astropy.org/Coordinates.html
and this:
https://docs.astropy.org/en/stable/coordinates/matchsep.html

and this function:
https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord.match_to_catalog_sky

"""

import astropy.table
import pandas as pd
import numpy as np

def match(tab, ds):
    """
    Matches a short input table `tab` with a long table/spreadsheet `ds`

    Optimal use case: you care a lot about seeing every possible match from `tab` into
    `ds`, and don't mind if the output is somewhat verbose.

    """
    # tab: input table with columns of ra, dec in decimal degrees (for now)
    # ds: "summary spreadsheet" pandas dataframe

    # simply print the separation between the two

    # for each item in tab:
    for item in tab:

        #    compute a distance between `item` and each element of `ds`.

        # for a locally flat approximation, the "rectangular" distance between two items 
        # in RA, Dec space has a term `delta` which is the cosine of the declination.

        # one degree longitudinal separation on the equator is longer than 
        #   one degree longitudinal separation near the poles.

        delta = np.cos(item.dec)

        dist = delta * (item.ra - ds.ra)**2 + ()



# c = SkyCoord(...)
# for i in range(len(c)):
#     ra, dec = c.ra[i], c.dec[i]
#     # do stuff

if __name__ == "__main__":

    # Here's a dummy 'test' example.

    tab = astropy.table.Table()
    tab_ra = [1, 2, 3]
    tab_dec = [1, 2, 3]

    tab.add_columns([tab_ra, tab_dec], names=["RA", "Dec"])

    ds = pd.DataFrame()
    ds.    
