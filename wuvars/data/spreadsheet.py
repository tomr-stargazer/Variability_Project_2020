"""
This is a place to import data from. Here's where we can grab the summary spreadsheets, 
as data structures.

"""

import os
from datetime import datetime
from astropy.table import Table
import pandas as pd

from recordclass import recordclass

wserv_ids = [1, 5, 7, 8, 11]
wserv_ids_v2 = [5, 7, 8, 11]

Dataset = recordclass(
    "Dataset", (f"wserv{n}" for n in wserv_ids), defaults=(None,) * len(wserv_ids)
)

Dataset_v2 = recordclass(
    "Dataset", (f"wserv{n}" for n in wserv_ids_v2), defaults=(None,) * len(wserv_ids_v2)
)

def load_v1():

    v1 = Dataset()

    spreadsheet_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/analysis_artifacts"

    spreadsheets = []

    for i, wserv in enumerate(wserv_ids):
        spreadsheet_path = os.path.join(
            spreadsheet_root,
            f"wserv{str(wserv)}",
            f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_new_error_corrected_summary_spreadsheet.h5",
        )

        ds = pd.read_hdf(spreadsheet_path, key="table")
        # adjust post-facto for the missing 2/3 weight.
        ds.loc[:, ("variability", "Stetson_JHK")] *= 1 / 1.5

        v1[i] = ds

    return v1


def load_wserv_v2(wserv):
    """
    Loads a photometry data given a WSERV id.

    """
    spreadsheet_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/analysis_artifacts"
    spreadsheet_path = os.path.join(
        spreadsheet_root,
        f"wserv{str(wserv)}",
        f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_new_error_corrected_summary_spreadsheet.h5",
    )

    ds = pd.read_hdf(spreadsheet_path, key="table")
    # adjust post-facto for the missing 2/3 weight.
    ds.loc[:, ("variability", "Stetson_JHK")] *= 1 / 1.5

    return ds


def load_v2():
    v2 = Dataset_v2()
    for i, wserv in enumerate(wserv_ids_v2):
        v2[i] = load_wserv_v2(wserv)

    return v2
