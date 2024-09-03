"""
Converts Spectral Type to effective Temperature, and vice versa, for young M-L dwarfs.

Explicitly uses the formulations of Luhman et al. (2003)
https://ui.adsabs.harvard.edu/abs/2003ApJ...593.1093L/abstract

and Faherty et al. (2016),
https://ui.adsabs.harvard.edu/abs/2016ApJS..225...10F/abstract
as compiled by Robberto et al. (2020).
https://ui.adsabs.harvard.edu/abs/2020ApJ...896...79R/abstract

Spectral type numbers start at M0 = 0; see `spectral_type_to_number.py` in this repo.

INPUT TABLES
------------

Luhman et al. (2003):
TABLE 8
Adopted Temperature Scale
Spectral Type ............... Teff (K)

M1 .............................. 3705
M2 .............................. 3560
M3 .............................. 3415
M4 .............................. 3270 
M5 .............................. 3125 
M6 .............................. 2990 
M7 .............................. 2880 
M8 .............................. 2710 
M9 .............................. 2400

Robberto et al. (2020)
Table A2
Spectral Type vs. Effective Temperature Relation 
According to Faherty et al. (2016)

Spectral Type ..... Teff(K)

L0 ................... 2173
L1 ................... 1983
L2 ................... 1809
L3 ................... 1649
L4 ................... 1504

"""

import numpy as np
from wuvars.analysis.spectral_type_to_number import get_num_from_SpT, get_SpT_from_num

# Luhman 2003
spt_teff_pairs = [
    ("M1", 3705),
    ("M2", 3560),
    ("M3", 3415),
    ("M4", 3270),
    ("M5", 3125),
    ("M6", 2990),
    ("M7", 2880),
    ("M8", 2710),
    ("M9", 2400),
]

# Herzceg & Hillenbrand 2014
hh14_spt_teff_pairs = [
    ("M1", 3720),
    ("M2", 3560),
    ("M3", 3410),
    ("M4", 3190),
    ("M5", 2980),
    ("M6", 2860),
    ("M7", 2770),
    ("M8", 2670),
    ("M9", 2570),
]

# Faherty 2016
Faherty_spt_teff_pairs = [
    ("L0", 2173),
    ("L1", 1983),
    ("L2", 1809),
    ("L3", 1649),
    ("L4", 1504),
]

all_spt_teff_pairs = spt_teff_pairs + Faherty_spt_teff_pairs

all_spts = []
all_teffs = []
for spt, teff in all_spt_teff_pairs:
    all_spts.append(spt)
    all_teffs.append(teff)

spt_num_array = np.array([get_num_from_SpT(x) for x in all_spts])
teff_array = np.array(all_teffs)


hh14_all_spt_teff_pairs = hh14_spt_teff_pairs + Faherty_spt_teff_pairs

hh14_all_spts = []
hh14_all_teffs = []
for spt, teff in hh14_all_spt_teff_pairs:
    hh14_all_spts.append(spt)
    hh14_all_teffs.append(teff)

hh14_spt_num_array = np.array([get_num_from_SpT(x) for x in hh14_all_spts])
hh14_teff_array = np.array(hh14_all_teffs)


def get_Teff_from_SpT(SpT):
    """ 
    Uses simple interpolation to compute Teff for any SpT. 

    Takes strings `"M8"` or numeric values `8`. "M0" is 0; "L0" is 10.

    """

    if isinstance(SpT, (str, np.str_)):
        SpT_num = get_num_from_SpT(SpT)
    else:
        SpT_num = SpT

    teff = np.interp(SpT_num, spt_num_array, teff_array)

    return teff


def get_SpT_from_Teff(Teff, out="num"):

    SpT_num = np.interp(Teff, teff_array[::-1], spt_num_array[::-1])

    if out == "num":
        return SpT_num
    else:
        SpT = get_SpT_from_num(SpT_num)
        return SpT


def get_Teff_from_SpT_HH14(SpT):
    """ 
    Uses simple interpolation to compute Teff for any SpT. 

    Takes strings `"M8"` or numeric values `8`. "M0" is 0; "L0" is 10.

    """

    if isinstance(SpT, (str, np.str_)):
        SpT_num = get_num_from_SpT(SpT)
    else:
        SpT_num = SpT

    teff = np.interp(SpT_num, hh14_spt_num_array, hh14_teff_array)

    return teff


def get_SpT_from_Teff_HH14(Teff, out="num"):

    SpT_num = np.interp(Teff, hh14_teff_array[::-1], hh14_spt_num_array[::-1])

    if out == "num":
        return SpT_num
    else:
        SpT = get_SpT_from_num(SpT_num)
        return SpT


# poor-man's tests.
if __name__ == "__main__":

    np.testing.assert_equal(get_Teff_from_SpT("M1"), 3705)
    np.testing.assert_equal(get_Teff_from_SpT(np.str_("M1")), 3705)
    np.testing.assert_equal(get_Teff_from_SpT(1), 3705)

    np.testing.assert_equal(get_SpT_from_Teff(2400, out="str"), "M9.0")
    np.testing.assert_equal(get_SpT_from_Teff(2400, out="num"), 9)

    # assert get_num_from_SpT("M0") == 0
    # assert get_num_from_SpT("M5.5") == 5.5
    # assert get_num_from_SpT("L1") == 11
    # assert get_num_from_SpT("K9") == -1
    # assert get_num_from_SpT("G2") == -18

    # assert get_SpT_from_num(0) == "M0"
    # assert get_SpT_from_num(5.5) == "M5.5"
    # assert get_SpT_from_num(11) == "L1"
    # assert get_SpT_from_num(-1) == "K9"
    # assert get_SpT_from_num(-18) == "G2"
