"""
A function to revise errorbars.

There are two conceivable approaches:
- a single function that takes a filename, does the thing, writes a file, returns nothing
- two functions:
-   a function that takes just an astropy table and either alters it in-place or returns a revised one
-   a function that wraps the above fn with file read and write

The second approach is more testable. But this is SO simple.

"""

import os
import pathlib
from datetime import datetime

import numpy as np
from astropy.table import Table
from wuvars.analysis.rederive_error_correction import decorrect_error


def revise_errorbars(improperly_error_corrected_dataset, s, c):

    dat = improperly_error_corrected_dataset

    bands = ["J", "H", "K"]
    for b in bands:

        # undo the old error correction terms

        improper_err = dat[f"{b}APERMAG3ERR"]
        err = decorrect_error(improper_err)

        # apply the new error correction terms
        new_corrected_err = np.sqrt(c * err ** 2 + s ** 2)

        dat[f"{b}APERMAG3ERR"] = new_corrected_err

    return dat


def load_and_revise_errorbars_and_save(input_path, output_path, s, c):
    print(f"Reading data from {input_path}")
    dat = Table.read(input_path)
    revised_dat = revise_errorbars(dat, s, c)
    print(f"Writing data to {output_path}")
    revised_dat.write(output_path, overwrite=True)

    return None


def load_and_revise_errorbars_and_save_debugger(input_path, output_path, s, c, rows=10):
    print(f"Reading data from {input_path}")
    dat = Table.read(input_path)[:10]
    revised_dat = revise_errorbars(dat, s, c)
    print(f"Writing data to {output_path}")
    revised_dat.write(output_path, overwrite=True)

    return None


def small_table_for_error_reproduction():
    # do this for each cleaned dataset.

    wserv_ids = [11]

    # Derived in the following script:
    # /Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/Fitting the error-correction terms to each dataset.ipynb
    error_correction_terms_s = [0.0093]
    error_correction_terms_c = [0.859]

    # where do these data live?
    input_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"
    # output_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/analysis_artifacts"

    for wserv, s, c in zip(
        wserv_ids, error_correction_terms_s, error_correction_terms_c
    ):
        input_path = os.path.join(
            input_root,
            f"wserv{str(wserv)}",
            f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5.h5",
        )

        output_path = os.path.join(
            input_root,
            f"wserv{str(wserv)}",
            f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_new_error_corrected_debug_small.h5",
        )

        pathlib.Path(os.path.join(input_root, f"wserv{str(wserv)}")).mkdir(
            parents=True, exist_ok=True
        )

        print(f"INPUT / OUTPUT for WSERV{wserv}:", input_path, output_path)
        startTime = datetime.now()
        print(f"Starting at: {startTime}")

        load_and_revise_errorbars_and_save_debugger(input_path, output_path, s, c)

        print(f"WSERV{wserv} elapsed time: ", datetime.now() - startTime)


if __name__ == "__main__":

    # do this for each cleaned dataset.

    wserv_ids = [1, 5, 7, 8, 11]

    # Derived in the following script:
    # /Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/Fitting the error-correction terms to each dataset.ipynb
    error_correction_terms_s = [0.0082, 0.0101, 0.0082, 0.0085, 0.0093]
    error_correction_terms_c = [0.921, 0.932, 0.919, 0.921, 0.859]

    # where do these data live?
    input_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"
    # output_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/analysis_artifacts"

    for wserv, s, c in zip(
        wserv_ids, error_correction_terms_s, error_correction_terms_c
    ):
        input_path = os.path.join(
            input_root,
            f"wserv{str(wserv)}",
            f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5.h5",
        )

        output_path = os.path.join(
            input_root,
            f"wserv{str(wserv)}",
            f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_new_error_corrected.fits",
        )

        pathlib.Path(os.path.join(input_root, f"wserv{str(wserv)}")).mkdir(
            parents=True, exist_ok=True
        )

        print(f"INPUT / OUTPUT for WSERV{wserv}:", input_path, output_path)
        startTime = datetime.now()
        print(f"Starting at: {startTime}")

        load_and_revise_errorbars_and_save(input_path, output_path, s, c)

        print(f"WSERV{wserv} elapsed time: ", datetime.now() - startTime)
