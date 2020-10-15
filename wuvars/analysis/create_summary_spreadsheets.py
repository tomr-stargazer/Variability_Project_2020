# For each graded_cleaned_scrubbed_clipped_whatever...
# run variability_selection.spreadsheet_maker

import os
import pathlib
from wuvars.analysis.variability_selection import spreadsheet_maker


def write_summary_spreadsheet(filename, output):
    # takes in a dataset from a given filename
    # makes an intermediate summary spreadsheet
    # writes that to file
    # returns nothing

    # load in the .fits.gz file - either via astropy or directly into pandas

    pass


if __name__ == "__main__":

    wserv_ids = [1, 5, 7, 8, 11]

    input_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"
    output_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/analysis_artifacts"

    for wserv in wserv_ids:
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

        pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)

        print(input_path, output_path)
        write_summary_spreadsheet(input_path, output_path)
