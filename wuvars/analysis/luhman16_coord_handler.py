"""
A quick little module to extract coordinates from the Luhhman et al. 2016 tables,
which provide separate columns for each sexagesimal component of the coordinates.

"""

import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord


def coords_from_Luhman_table(table):

    # RA strings like "HH MM SS.S..." with native precision
    ra_str = [
        f"{int(h):02d} {int(m):02d} {s}"
        for h, m, s in zip(table["RAh"], table["RAm"], table["RAs"])
    ]

    # Dec strings like "+DD MM SS.S..." with native precision
    dec_str = [
        f"{sign}{int(d):02d} {int(m):02d} {s}"
        for sign, d, m, s in zip(table["DE-"], table["DEd"], table["DEm"], table["DEs"])
    ]

    # Use SkyCoord directly on the sexagesimal strings
    coordinates_array = SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg))

    return coordinates_array
