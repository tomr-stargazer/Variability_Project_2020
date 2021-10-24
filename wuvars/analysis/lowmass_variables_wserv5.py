"""
This script has the formal definitions of all the low-mass variables in WSERV5 (ONC).

These variables consist of the following objects:

- low-mass stars (defined, here, by being in Robberto+20 and having a Teff<3200)
- which were approved by a by-eye inspection (see `bd_matching_v2`), and are either:
- periodic variables, or
- non-periodic variables.

Periodic variables were identified by a process

"""

import os
import numpy as np
import pandas as pd

# First, load up the periodic variables for WSERV5.

spreadsheet_dir = (
    "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes"
)
flags = ["Y", "Yw", "N", "YfY", "?fY", "YfYw", "?fYw", "YfN", "?fN"]
periodic_flags = [flag for flag in flags if flag[-1] in ("Y", "w")]
# print(periodic_flags)


def load_periodic_variables_wserv5(verbose=True):

    period_sheet_onc = pd.read_excel(
        os.path.join(spreadsheet_dir, "ONC_source_properties_periods_inspected.xlsx")
    )

    periodic_onc = period_sheet_onc[
        np.in1d(period_sheet_onc["Periodic?"], periodic_flags)
    ]

    if verbose:
        print(f"Total objects in ONC: {len(period_sheet_onc)}")

        n_periodic = np.sum(np.in1d(period_sheet_onc["Periodic?"], periodic_flags))
        print(f"Total periodic objects: {n_periodic}")
        print(f"Fraction periodic: {n_periodic/len(period_sheet_onc) * 100:.1f}%")

    # is this what we want to return?
    return period_sheet_onc, periodic_onc


def load_nonperiodic_variables_wserv5(verbose=True):
    pass


def load_all_variables_wserv5(verbose=True):

    # get all of the objects

    # get the periodic variables

    # get the nonperiodic variables

    # get the union of the two sets

    # return the union of those two sets

    pass


if __name__ == "__main__":
    load_periodic_variables_wserv5()
    load_nonperiodic_variables_wserv5()
