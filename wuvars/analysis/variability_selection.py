"""
This is a script that includes some notes on variability selection.

I'm hoping to nail down some semantics around how variables *should* be selected. 
There will be some troubleshooting here.

"""

from collections import namedtuple

import numpy as np

# one idea:
# given a pandas dataframe that's indexed by SOURCEID, MEANMJDOBS etc,
# do a bunch of systematic groupby / aggregates of it.

# A note: we'll want to pay close attention to how the sentinel values are handled
# (what values do they take? do we want to cast them to NaN right out the gate?
#  will that be handled here in the stats calculations or will the NaN-type values be
#  changed in-place at an earlier stage of the data processing? etc)

# Tom, you need to do this! (Okay, I think we figured out an approach here that works for me.)
# We're gonna do laundry now. [ ] and maybe a walk?

Spreadsheet = namedtuple("mean", "median", "std", "min", "max", "range")


def spreadsheet_maker(df):
    """
    Compute several summary properties of EVERY column in this table, grouped by SOURCEID:
    - mean
    - median
    - std / rms
    - min
    - max
    - range

    and also the following columns:
    - J_counts
    - H_counts
    - K_counts
    (possibly sorted by "good/noflag", "warn", "error" / from the ppErrBits)
    (see, for example, how )

    and some specifically variability-related columns:
    - chisq (J, H, K)
    - Stetson (JHK, JK, JH, HK)

    possibly with "robust" / "outlier-removed" versions of some or all of these.

    I'd like to perhaps group all these into a hierarchical data structure that looks like this:

    # ss : stands for SpreadSheet 
    spreadsheet.mean['JAPERMAG3']
    spread


    """

    # let's prototype this!
    # what does it look like to do this for a single data array?

    df.groupby("SOURCEID").apply(np.nanmean)

    summary_functions = [np.nanmean, np.nanmedian]

    return this_goofy_namedtuple


# Here's another thing we want: the ability to select variable stars, given the above "spreadsheet" / summary properties.
# Basically, I want to be able to tell my data:
# Q2 variables
# ------------
# - N_J, N_H, N_K > 50 (each)
# AND
# - JPPERRMAX, HPPERRMAX, KPPERRMAX = 0 (each)
# AND (any of the following)
# [
# - JHK Stetson > 1 OR
# - JK, HK, JH Stetson > 1 OR
# - chisq_J > N_J   OR
# - chisq_H > N_H   OR
# - chisq_K > N_K
# ]

# similar (but ultimately more convoluted) criteria can be devised for what we might call "Q1" variables, which
# would be actually broken into "two-band" Q1 variables, and "one-band" Q1 variables:
# - JH, HK, JK
# - J, H, K
# more can be written soon; a lot of it depends on how many

# There's a slight philosophical distinction between the above framework and how I think the WSA does things.
# In my framework, I disqualify a star from a quality grade if it has a *single* flagged night


def select_variables(spreadsheet, parameters):
    """
    Selects variables from a spreadsheet, given parameters.

    Parameters:
    - minimum Stetson that qualifies variability
    - 

    Returns: ???
    - Possibly an index array that maps back to the input spreadsheet? like a boolean array
      so that you could say 
      `variable_sourceids = spreadsheet['SOURCEIDS'][variable_indices]`
      given the variable_indices return value

    Notes:
    - I'd like to make it straightforward to ... do unions and intersections of different variable subsets.

    """
