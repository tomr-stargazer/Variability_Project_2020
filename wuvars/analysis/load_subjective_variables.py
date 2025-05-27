"""
This is a module that you import in order to access the subjective variables (V0).

It borrows heavily from `load_periodics_v4.py`.

"""

import os

import numpy as np
import pandas as pd
from wuvars.data import spreadsheet

subjective_results_dir = (
    "/Users/tsrice/Documents/Variability_Project_2020/Results/Q0_followup_v4/"
)

ngc_subjectives = pd.read_excel(
    os.path.join(
        subjective_results_dir,
        "ngc",
        "Q0_objects_for_followup_ngc_2024_Oct_13_inspected.xlsx",
    )
)

ic_subjectives = pd.read_excel(
    os.path.join(
        subjective_results_dir,
        "ic",
        "Q0_objects_for_followup_ic_2024_Oct_13_inspected.xlsx",
    )
)

subjectives_dict = {7: ngc_subjectives, 8: ic_subjectives}

def select_subjective_variables(wserv):

    subjectives_sheet = subjectives_dict[wserv]

    subjectives = subjectives_sheet[subjectives_sheet['VARIABLE'] == 1]

    # periodics = period_sheet[np.in1d(period_sheet["Periodic?"], periodic_flags)]
    subjectives_by_sourceid = subjectives.set_index("SOURCEID")

    return subjectives_by_sourceid
