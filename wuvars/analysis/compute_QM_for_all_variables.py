"""
We're going to compute and save the Q, M values (in each of the three bands)
for all variables, and save the output so that we don't have to re-compute.

"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wuvars.analysis.variability_selection as sv
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.calculate_M import compute_M
from wuvars.analysis.calculate_Q import compute_Q_automatically
from wuvars.analysis.spectral_type_to_number import get_SpT_from_num
from wuvars.analysis.variability_selection_curved import (curve_Stetson, sv_hk,
                                                          sv_jh, sv_jhk, sv_jk)
from wuvars.data import photometry, quality_classes, spreadsheet

warnings.filterwarnings("ignore")

# This part pulls up the source information from our position-matching scheme.
ngc_match = match_ngc()
ic_match = match_ic()

names = ["ngc", "ic"]
wserv_dict = {"ngc": 7, "ic": 8}
fullname_dict = {"ngc": "NGC 1333", "ic": "IC 348"}
match_dict = {"ngc": ngc_match, "ic": ic_match}
spread_dict = {"ngc": spreadsheet.load_wserv_v2(7), "ic": spreadsheet.load_wserv_v2(8)}
q_dict = {"ngc": quality_classes.load_q(7), "ic": quality_classes.load_q(8)}
photometry_dict = {
    "ngc": photometry.group_wserv_v2(
        photometry.load_wserv_v2(7, suffix="_outlier_cleaned_152")
    ),
    "ic": photometry.group_wserv_v2(photometry.load_wserv_v2(8)),
}


# load up all variables
# which is necessarily going to be split into the two regions...
# for region in regions:
#   for variable in variables:

colnames = [
    "SOURCEID",
    "M_J",
    "M_H",
    "M_K",
    "Q_J",
    "Per_J",
    "Q_H",
    "Per_H",
    "Q_K",
    "Per_K",
]


if __name__ == "__main__":
    # Let's do the overview of IC 348 and NGC 1333.
    # How many objects were in the original catalog
    print("The input catalogs (total):\n")
    for name in names:
        print("***")
        print(fullname_dict[name], "\n")
        match = match_dict[name]
        spread = spread_dict[name]
        dataset = photometry_dict[name]

        region_data = []

        for i, sid in enumerate(match.approved["SOURCEID"]):

            M_dict = compute_M(dataset, sid)
            Q_dict = compute_Q_automatically(dataset, sid, verbose=False)

            # print(sid)
            # print(M_dict)
            # print(Q_dict)

            this_row = [
                sid,
                M_dict["J"],
                M_dict["H"],
                M_dict["K"],
                *Q_dict["J"],
                *Q_dict["H"],
                *Q_dict["K"],
            ]
            region_data.append(this_row)

            # print(this_row)

        region_df = pd.DataFrame(region_data, columns=colnames)
        region_df.to_hdf(f"{name}_QM.h5", key=f"{name}_QM")
        print(f"Wrote `{name}_QM.h5` to file!")


# dataset will look like
# dat11 = photometry.group_wserv_v2(photometry.load_wserv_v2(11))
