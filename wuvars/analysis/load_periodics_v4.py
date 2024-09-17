"""
This is a module that you import in order to access the v4 periods.

"""

import os

import numpy as np
import pandas as pd
from wuvars.data import spreadsheet

prototypes_dir = (
    "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes"
)

ngc_periods = pd.read_excel(
    os.path.join(prototypes_dir, "v4_NGC_source_properties_periods_inspected.xlsx")
)

ic_periods = pd.read_excel(
    os.path.join(prototypes_dir, "v4_IC_source_properties_periods_inspected.xlsx")
)

periods_dict = {7: ngc_periods, 8: ic_periods}


for _periods in [ngc_periods, ic_periods]:

    for i, row in _periods.iterrows():
        if row["Periodic?"] == "Y":
            band = row["Best Band"]
            method = row["Method"]

            per_col = f"{method}_period_{band}"
            amp_col = f"{method}_per_amp_{band}"

            # loc[:, ('Period', i)]
            _periods['Period'][i] = row[per_col]
            _periods["Amp"][i] = row[amp_col]

# now ngc_periods and ic_periods are importable.

def select_periodic_variables_v4(wserv):

    period_sheet = periods_dict[wserv]

    periodics = period_sheet[period_sheet['Periodic?'] == 'Y']

    # periodics = period_sheet[np.in1d(period_sheet["Periodic?"], periodic_flags)]
    periodics_by_sourceid = periodics.set_index("SOURCEID")

    return periodics_by_sourceid
