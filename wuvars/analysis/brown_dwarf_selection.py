"""
Here's where I'll try to write code that selects brown dwarfs based on their colors.

"""
from wuvars.analysis.bd_mags import apparent_BD_mags_jhk

import pdb


def simple_BD_select(ds, wserv):

    # rieke & lebofsky 1985 reddening law
    A_K = 0.112
    A_H = 0.175
    E_HK = A_H - A_K
    slope = A_K / E_HK  # this is like a slope in color-mag space

    J_BDlimit = apparent_BD_mags_jhk[wserv][0]
    H_BDlimit = apparent_BD_mags_jhk[wserv][1]
    K_BDlimit = apparent_BD_mags_jhk[wserv][2]

    HK_BDlimit = H_BDlimit - K_BDlimit

    K_intercept = K_BDlimit - slope * HK_BDlimit

    h = ds["median"]["HAPERMAG3"]
    k = ds["median"]["KAPERMAG3"]
    hmk = h - k

    bd_selection = (hmk > HK_BDlimit) & (k > slope * hmk + K_intercept)

    return bd_selection


def simple_BD_select_v2(ds, wserv):

    # rieke & lebofsky 1985 reddening law
    A_K = 0.112
    A_H = 0.175
    E_HK = A_H - A_K
    slope = A_K / E_HK  # this is like a slope in color-mag space

    J_BDlimit = apparent_BD_mags_jhk[wserv][0]
    H_BDlimit = apparent_BD_mags_jhk[wserv][1]
    K_BDlimit = apparent_BD_mags_jhk[wserv][2]

    HK_BDlimit = H_BDlimit - K_BDlimit

    K_intercept = K_BDlimit - slope * HK_BDlimit

    h = ds["median"]["HAPERMAG3"]
    k = ds["median"]["KAPERMAG3"]
    hmk = h - k

    bd_selection = (hmk > HK_BDlimit) & (k > slope * hmk + K_intercept) | (
        hmk <= HK_BDlimit
    ) & (k >= K_BDlimit)

    return bd_selection


# brainstorming something more complex:
"""

def less_simple_BD_select(ds, wserv):

    # rieke & lebofsky 1985 reddening law
    A_K = 0.112
    A_H = 0.175
    E_HK = A_H - A_K  
    slope = A_K / E_HK # this is like a slope in color-mag space

    J_BDlimit = apparent_BD_mags_jhk[wserv][0]
    H_BDlimit = apparent_BD_mags_jhk[wserv][1]
    K_BDlimit = apparent_BD_mags_jhk[wserv][2]

    HK_BDlimit = H_BDlimit - K_BDlimit

    K_intercept = K_BDlimit - slope * HK_BDlimit

    j = ds["median"]["JAPERMAG3"]
    h = ds["median"]["HAPERMAG3"]
    k = ds["median"]["KAPERMAG3"]
    hmk = h - k


    bd_selection_color = (hmk > HK_BDlimit) & (k > slope * hmk + K_intercept)

    # also include stars with NAN as their H and J mags 

    bd_selection_limit = (np.isnan(h) & np.isnan(j) & ~np.isnan(k))

    bd_selection = bd_selection_color | bd_selection_limit

    return bd_selection
"""
