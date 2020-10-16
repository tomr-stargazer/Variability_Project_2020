# For each graded_cleaned_scrubbed_clipped_whatever...
# run variability_selection.spreadsheet_maker

import os
from datetime import datetime
import pathlib
from astropy.table import Table

from wuvars.analysis.variability_selection import spreadsheet_maker


def write_summary_spreadsheet(filename, output):
    # takes in a dataset from a given filename
    # makes an intermediate summary spreadsheet
    # writes that to file
    # returns nothing

    # load in the .fits.gz (or hdf5) file - via astropy, then export to pandas

    dat = Table.read(filename)
    df = dat.to_pandas()

    ds = spreadsheet_maker(df)
    ds.to_hdf(output, key="table")

    return None


if __name__ == "__main__":

    wserv_ids = [1, 5, 7, 8, 11]

    input_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"
    output_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/analysis_artifacts"

    for wserv in wserv_ids[::-1]:
        input_path = os.path.join(
            input_root,
            f"wserv{str(wserv)}",
            f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5.h5",
        )

        output_path = os.path.join(
            output_root,
            f"wserv{str(wserv)}",
            f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_summary_spreadsheet.h5",
        )

        pathlib.Path(os.path.join(output_root, f"wserv{str(wserv)}")).mkdir(
            parents=True, exist_ok=True
        )

        print(f"INPUT / OUTPUT for WSERV{wserv}:", input_path, output_path)
        startTime = datetime.now()
        print(f"Starting at: {startTime}")

        write_summary_spreadsheet(input_path, output_path)

        print(f"WSERV{wserv} elapsed time: ", datetime.now() - startTime)
