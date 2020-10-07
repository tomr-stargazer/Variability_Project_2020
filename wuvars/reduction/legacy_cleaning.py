"""
This script handles data cleaning, using legacy code.

I'm not proud of that legacy code, but - so it goes.

Depends on ATpy, because of (aforementioned) legacy reasons.

Also depends on all code from the following projects being importable
(as individual modules! not even as submodules within projects!) on the pythonpath:
- https://github.com/tomr-stargazer/wuvars-proto
- https://github.com/tomr-stargazer/wuvars-orion

It's not up to par with modern python packaging standards (tbh, it never was). Needs
refactoring badly. But it does... "technically" work, when set up properly.

"""

import os
import numpy as np
import atpy

from wuvars.reduction.legacy_cleaning_helpers import spreadsheet_mapreduce_combine
import night_cleanser
import variability_map

# let's pseudocode out the steps, with perseus_data_processing as a partial guide

# load up the .fits file
# fix the error bars
# clean the data:
# - first round of spreadsheet-making as input to the data cleaning (finding constants)
#     input: the error-corrected data
#     output: summary spreadsheet
# - variability_map.exposure_grader  - makes these jvc things
#     input: summary spreadsheet
#     output: hash of timestamps -> quality grades
# - use the quality grades to cleanse, scrub, and dust
#   (actually three separate steps)
#   - cleanse with null_cleanser_grader
#     input: photometry data, quality grades, a minimum threshold
#     output: "cleansed" photometry data
#   - scrub with selective_flag_scrubber
#     input: cleansed_data, "minimum" spreadsheet
#     output: cleansed_scrubbed_data
#   - dust with errorbar_duster
#     input: cleansed_scrubbed data
#     output: cleansed_scrubbed_dusted data
# - save the output after all that!
#
# we want a "makefile" or something similar that builds all of the above, ideally with
# a single, repeatable command.


def do_it(
    input_location,
    output_location,
    name,
    threshold=0.95,
    clobber=True,
    minimum_nights=50,
    splits=10,
):
    """
    The idea here is that we get an input, we put the end product in the output.

    """

    print("Input: ", input_location)
    print("Output: ", output_location)
    print("Name: ", name)

    # load up the .fits file
    data = atpy.Table(input_location)
    data.table_name = name

    # fix the error bars
    s = 0.021
    c = 1.082

    data.JAPERMAG3ERR = np.sqrt(c * data.JAPERMAG3ERR ** 2 + s ** 2)
    data.HAPERMAG3ERR = np.sqrt(c * data.HAPERMAG3ERR ** 2 + s ** 2)
    data.KAPERMAG3ERR = np.sqrt(c * data.KAPERMAG3ERR ** 2 + s ** 2)

    data.JMHPNTERR = np.sqrt(c * data.JMHPNTERR ** 2 + s ** 2)
    data.HMKPNTERR = np.sqrt(c * data.HMKPNTERR ** 2 + s ** 2)

    # clean the data:
    # - first round of spreadsheet-making as input to the data cleaning (finding constants)
    #     input: the error-corrected data
    #     output: summary spreadsheet

    intermediate_location = os.path.join(output_location, "intermediate")
    from pathlib import Path

    Path(intermediate_location).mkdir(parents=True, exist_ok=True)

    try:
        summary_spreadsheet = atpy.Table(
            os.path.join(output_location, "summary_spreadsheet.h5")
        )
    except FileNotFoundError:
        summary_spreadsheet = spreadsheet_mapreduce_combine(data, intermediate_location, splits=splits)
        summary_spreadsheet.write(
            os.path.join(output_location, "summary_spreadsheet.h5"), overwrite=True
        )

    # spreadsheet = atpy.Table(path+"DATA/spreadsheet/full_data_errorcorrected_ce_spreadsheet.fits")
    minimum = summary_spreadsheet.where(
        (summary_spreadsheet.N_j >= minimum_nights)
        | (summary_spreadsheet.N_k >= minimum_nights)
        | (summary_spreadsheet.N_h >= minimum_nights)
    )
    minimum_constant = minimum.where(
        (minimum.j_rmsr < 0.05) & (minimum.h_rmsr < 0.05) & (minimum.k_rmsr < 0.05)
    )
    constants = minimum_constant.where(
        minimum_constant.Stetson < np.median(minimum_constant.Stetson)
    )
    # - variability_map.exposure_grader  - makes these jvc things

    #     input: summary spreadsheet
    #     output: hash of timestamps -> quality grades
    try:
        jvc_ratio = np.loadtxt(os.path.join(intermediate_location, "jvc_ratio.txt"))
        hvc_ratio = np.loadtxt(os.path.join(intermediate_location, "hvc_ratio.txt"))
        kvc_ratio = np.loadtxt(os.path.join(intermediate_location, "kvc_ratio.txt"))
    except OSError:
        jvc = variability_map.exposure_grader(data, constants, "j", 17)
        hvc = variability_map.exposure_grader(data, constants, "h", 16.7)
        kvc = variability_map.exposure_grader(data, constants, "k", 16)

        jvc_ratio = jvc[2]
        hvc_ratio = hvc[2]
        kvc_ratio = kvc[2]

        np.savetxt(os.path.join(intermediate_location, "jvc_ratio.txt"), jvc_ratio)
        np.savetxt(os.path.join(intermediate_location, "hvc_ratio.txt"), hvc_ratio)
        np.savetxt(os.path.join(intermediate_location, "kvc_ratio.txt"), kvc_ratio)

    # - use the quality grades to cleanse, scrub, and dust
    #   (actually three separate steps)
    #   - cleanse with null_cleanser_grader
    #     input: photometry data, quality grades, a minimum threshold
    #     output: "cleansed" photometry data
    #   - scrub with selective_flag_scrubber
    #     input: cleansed_data, "minimum" spreadsheet
    #     output: cleansed_scrubbed_data
    #   - dust with errorbar_duster
    #     input: cleansed_scrubbed data
    #     output: cleansed_scrubbed_dusted data

    timestamps = np.sort(list(set(data.MEANMJDOBS)))

    new_data = atpy.Table()
    for col in data.columns:
        new_data.add_column(col, data[col])

    cleansed_data = night_cleanser.null_cleanser_grader(
        new_data, timestamps, jvc_ratio, hvc_ratio, kvc_ratio, threshold=threshold
    )
    cleansed_scrubbed_data = night_cleanser.selective_flag_scrubber(
        cleansed_data, minimum
    )
    cleansed_scrubbed_dusted_data = night_cleanser.errorbar_duster(
        cleansed_scrubbed_data
    )

    # - save the output after all that!
    datalist = [cleansed_data, cleansed_scrubbed_data, cleansed_scrubbed_dusted_data]
    pathlist = [
        f"{name}_graded_clipped0.95.h5",
        f"{name}_graded_clipped0.95_scrubbed0.1.h5",
        f"{name}_graded_clipped0.95_scrubbed0.1_dusted0.5.h5",
    ]

    for datafile, path in zip(datalist, pathlist):

        if os.path.isfile(os.path.join(output_location, path)) and clobber:
            os.remove(os.path.join(output_location, path))
        datafile.write(os.path.join(output_location, path))


def make():

    raw_data_path = (
        "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/Raw_Downloads"
    )
    # do it for wserv1
    wserv1_data_input = os.path.join(raw_data_path, "wserv1.fits.gz")
    wserv1_output_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/wserv1/"
    do_it(wserv1_data_input, wserv1_output_path, "WSERV1", splits=100)

    # do it for wserv5
    wserv5_data_input = os.path.join(raw_data_path, "wserv5.fits.gz")
    wserv5_output_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/wserv5/"
    do_it(wserv5_data_input, wserv5_output_path, "WSERV5", splits=25)

    # do it for wserv7
    wserv7_data_input = os.path.join(raw_data_path, "wserv7.fits.gz")
    wserv7_output_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/wserv7/"
    do_it(wserv7_data_input, wserv7_output_path, "WSERV7", splits=25)

    # do it for wserv8
    wserv8_data_input = os.path.join(raw_data_path, "wserv8.fits.gz")
    wserv8_output_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/wserv8/"
    do_it(wserv8_data_input, wserv8_output_path, "WSERV8", splits=25)

    # do it for wserv11
    wserv11_data_input = os.path.join(raw_data_path, "wserv11.fits.gz")
    wserv11_output_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/wserv11/"
    do_it(wserv11_data_input, wserv11_output_path, "WSERV11", splits=25)

    # now we should have reduced, cleaned photometry for all five regions.
    pass


def test_make():
    # do it for a narrow subset of wserv7.
    # like, a quarter of the dates, just between mags 13 and 16.
    # let's download that now.
    path_to_prototype = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/Raw_Downloads/prototype/"
    prototype_input_file = os.path.join(
        path_to_prototype,
        "wserv7_nodupes_timeslice.lt56200_magslice.jhk.btw.13.16.fits.gz",
    )

    prototype_output_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/prototypes/"

    do_it(prototype_input_file, prototype_output_path, "DUMMY", minimum_nights=15)
