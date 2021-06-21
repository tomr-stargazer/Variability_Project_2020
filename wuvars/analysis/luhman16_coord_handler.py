"""
A quick little module to extract coordinates from the Luhhman et al. 2016 tables,
which provide separate columns for each sexagesimal component of the coordinates.

"""

import numpy as np
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u


def coords_from_Luhman_table(table):

    # get the coordinates into a usable state
    RA = Angle(
        (table["RAh"].data, table["RAm"].data, table["RAs"].data), unit=u.hourangle
    )
    Dec_sign = np.where(table["DE-"] == "+", 1, -1)
    Dec = Angle((Dec_sign * table["DEd"], table["DEm"], table["DEs"],), unit=u.deg,)

    coordinates_array = SkyCoord(ra=RA, dec=Dec)

    return coordinates_array
