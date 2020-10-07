"""
All of this is stolen, like, verbatim from `variability_script_perseus.py` from wuvars-perseus.

I might hack it to pieces, though.

Original docstring:

This is a script (not a program!) that uses my functions to 
quantify variability for my stars. 

To edit the spreadsheet parameters, C-s "EEEEE"

"""

import os
import datetime

import atpy

import spread3 as sp


print("Hey, just a heads-up, this is an INTERACTIVE script.")
print(" You should call the following functions:")
print(" -test() # To make sure everything's working fine")
print("         # before wasting a lot of time.")
print(" -calculate_stuff() # To calculate stuff.")
print(" -glue_stuff() # To glue together the calculated stuff.")
print("               # Note, this one returns the spreadsheet.")
print("")
print("New feature: you can pass a number to calculate_stuff() and glue_stuff()")
print("(such as 25, 50, 100) as a manual control on how many chunks to split")
print("the data into. Make sure to use the same number for both functions!!")


def spreadsheet_mapreduce(data, intermediate_location, splits=10, start=0):
    """
    Runs the spreadsheet, first splitting it into `splits` 
    spreadsheets and then joining them. 

    """

    if type(splits) is not int or type(start) is not int:
        raise TypeError

    # We are going to split this into 10 smaller pieces through the magic of mod operations! woo.

    split_data = []

    for i in range(start, splits):
        data_i = data.where(data.SOURCEID % splits == i)

        split_data.append(data_i)

        lookup_i = sp.base_lookup(data_i)

        # The parameter "-1" is the season that tells data_cut not to make
        # any cuts on the data.
        sp.spreadsheet_write(
            data_i,
            lookup_i,
            -1,
            os.path.join(intermediate_location, "sp%d.h5" % i),
            flags=256,
            per=False,
            graded=False,
            rob=True,
            colorslope=True,
        )
        # EEEEE this is a flag to come and find this section of code

        try:
            now = datetime.datetime.strftime(
                datetime.datetime.now(), "%Y-%m-%d %H:%M:%S"
            )
        except:
            now = "sometime"
        print(("finished chunk %d at %s" % (i, now)))


def spreadsheet_combine(intermediate_location, splits=10, start=0):
    """ Read in the tables from earlier and glue them together """

    if type(splits) is not int:
        raise TypeError

    spread = atpy.Table(os.path.join(intermediate_location, "sp%d.h5" % start))

    for i in range(1 + start, splits):
        other_spread = atpy.Table(os.path.join(intermediate_location, "sp%d.h5" % i))
        spread.append(other_spread)

    return spread


def spreadsheet_mapreduce_combine(
    data, intermediate_location, splits=10, combine_start=0, glue_start=0
):
    spreadsheet_mapreduce(
        data, intermediate_location, splits=splits, start=combine_start,
    )
    return spreadsheet_combine(intermediate_location, splits=splits, start=glue_start)
