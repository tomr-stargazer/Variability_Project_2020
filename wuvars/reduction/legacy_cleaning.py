"""
This script handles data cleaning, using legacy code.

I'm not proud of that legacy code, but - so it goes.

"""

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

def do_it(input_location, output_location, name):
    """
    The idea here is that we get an input, we put the end product in the output.

    """

    print("Input: ", input_location)
    print("Output: ", output_location)
    print("Name: ", name)


def make():

    # do it for wserv1

    # do it for wserv5

    # do it for wserv7

    # do it for wserv8

    # do it for wserv11

    # now we should have reduced, cleaned photometry for all five regions.
    pass


def test_make():
    # do it for a narrow subset of wserv7.
    # like, a quarter of the dates, just between mags 13 and 16.
    # let's download that now.
    pass