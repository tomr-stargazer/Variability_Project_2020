"""
This script actually applies the outlier cleaning to NGC 1333.

Pro tip: I am only applying this (currently) to the objects in ngc_match.approved.

If I do another project (maybe on, like, the protostars? Or field variables) 
in NGC 1333 I'll have to do this process on those guys separately. 

But this is an awkward and intensive enough process that I don't want to just run it 
on all objects in the photometric data set.

2 September 2024 - Tom Rice

"""

import os
import pathlib
from datetime import datetime

import numpy as np
from astropy.table import Table
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.reduction.final_outlier_clean import (calculate_diffs,
                                                  identify_outliers)

null = np.double(-9.99999488e08)


def nullify_outliers(dataset, source_list):

    subdat = dataset[np.in1d(dataset["SOURCEID"], source_list)]
    cleaned_dataset = subdat.copy(copy_data=True)

    subdat_grouped = subdat.group_by("SOURCEID")

    for sid in source_list:

        diffs = calculate_diffs(subdat_grouped, sid)
        outliers = identify_outliers(
            subdat_grouped, sid, diffs, date_offset=56141, verbose=False
        )

        for (band, b_outliers) in outliers.items():

            for b_date in b_outliers:
                # nullify this one

                outlier_selection = cleaned_dataset["MEANMJDOBS"] == b_date
                cleaned_dataset[f"{band}APERMAG3"][outlier_selection] = null

    return cleaned_dataset


if __name__ == "__main__":

    startTime = datetime.now()
    print(f"Starting at: {startTime}")

    # we're just doing this for NGC 1333 (WSERV7), and only for 152 sources.
    ngc_match = match_ngc()
    ngc_dat = photometry.group_wserv_v2(photometry.load_wserv_v2(7))

    wserv = 7
    suffix = "_outlier_cleaned_152"

    output_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"
    filename = f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_new_error_corrected{suffix}.fits"
    full_output = os.path.join(output_path, f"wserv{str(wserv)}", filename)
    test_output = os.path.join(output_path, f"wserv{str(wserv)}", f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_new_error_corrected_152.fits")

    subdat = ngc_dat[np.in1d(ngc_dat["SOURCEID"], ngc_match.approved["SOURCEID"])]
    subdat.write(test_output, overwrite=True)

    if False:
        # testing something.
        cleaned_dataset = nullify_outliers(ngc_dat, [44508746107200])
    if True:

        cleaned_dataset = nullify_outliers(ngc_dat, ngc_match.approved["SOURCEID"])
        print(f"Writing to {full_output}")
        cleaned_dataset.write(full_output, overwrite=True)

        print(f"WSERV{wserv} elapsed time: ", datetime.now() - startTime)
