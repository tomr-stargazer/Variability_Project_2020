"""
This module allows you to import a canonical set of quality selections for each 
"""

import os

import pandas as pd
from wuvars.analysis.variability_selection import (sq0, sq1, sq1_h, sq1_j,
                                                   sq1_k, sq2)
from wuvars.data.custom_recordclass import recordclass
from wuvars.data.spreadsheet import load_wserv_v2

wserv_ids = [1, 5, 7, 8, 11]
n_obs_list = [130, 120, 171, 85, 110]

n_min_list = [60, 35, 80, 55, 65]
n_max_list = [100, 90, 160, 80, 100]

n_min_dict = {wserv: n_min for wserv, n_min in zip(wserv_ids, n_min_list)}
n_max_dict = {wserv: n_max for wserv, n_max in zip(wserv_ids, n_max_list)}

Qualityset = recordclass("Qualityset", ["q2", "q1", "q1_j", "q1_h", "q1_k", "q0"])


def load_q(wserv):

    ds = load_wserv_v2(wserv)

    q0 = sq0(ds, n_min_dict[wserv], n_max_dict[wserv])
    q1 = sq1(ds, n_min_dict[wserv], n_max_dict[wserv])
    q2 = sq2(ds, n_min_dict[wserv], n_max_dict[wserv])

    q1_j = sq1_j(ds, n_min_dict[wserv], n_max_dict[wserv])
    q1_h = sq1_h(ds, n_min_dict[wserv], n_max_dict[wserv])
    q1_k = sq1_k(ds, n_min_dict[wserv], n_max_dict[wserv])

    return Qualityset(q2, q1, q1_j, q1_h, q1_k, q0)
