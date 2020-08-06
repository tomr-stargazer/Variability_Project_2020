"""
This is a script that includes some notes on variability selection.

I'm hoping to nail down some semantics around how variables *should* be selected. 
There will be some troubleshooting here.

"""

import numpy as np
import pandas as pd

# one idea:
# given a pandas dataframe that's indexed by SOURCEID, MEANMJDOBS etc,
# do a bunch of systematic groupby / aggregates of it.

# A note: we'll want to pay close attention to how the sentinel values are handled
# (what values do they take? do we want to cast them to NaN right out the gate?
#  will that be handled here in the stats calculations or will the NaN-type values be
#  changed in-place at an earlier stage of the data processing? etc)

# Tom, you need to do this! (Okay, I think we figured out an approach here that works for me.)
# We're gonna do laundry now. [X] and maybe a walk?

# let's write up the variability functions - I think I implemented these recently in
# stetson_2020.py. (in Variability_2019)
# From c. 2012.


def delta(m, sigma_m, mean_m, n):
    """
    Normalized residual / "relative error" for one observation.
    Used in Stetson's J variability index.

    INPUTS:
        m: a single magnitude measurement in a certain band
        sigma_m: the uncertainty on that measurement
        mean_m: the mean magnitude in that band
        n: the number of observations in that band

    OUTPUTS:
        d: the "relative error"
    """

    d = np.sqrt(n / (n - 1)) * (m - mean_m) / sigma_m

    return d


# From c. 2012.
def S_threeband(j, sigma_j, h, sigma_h, k, sigma_k):
    """
    Computes the Stetson variability index for one star that has
    3 observations on each night. Uses Carpenter et al.'s notation.
    Simplified from the general expression assuming 3 observations every
      night.

    INPUTS:
        j: an array of J-band magnitudes
        sigma_j: an array of corresponding J-band uncertainties
        h: an array of H-band magnitudes
        sigma_h: an array of corresponding H-band uncertainties
        k: an array of K-band magnitudes
        sigma_k: an array of corresponding K-band uncertainties

    OUTPUTS:
        s: the Stetson variability index for 3 bands

    """

    n = j.size

    # Perhaps hackish
    if n < 2:
        return 0

    d_j = delta(j, sigma_j, np.nanmean(j), n)
    d_h = delta(h, sigma_h, np.nanmean(h), n)
    d_k = delta(k, sigma_k, np.nanmean(k), n)

    P_i = np.array([d_j * d_h, d_h * d_k, d_j * d_k])

    # I originally had two sums going: one over P_i, and one over all the
    # elements of n, but then I realized that a single sum over all axes
    # did the exact same thing (I tested it) so now it's back to one sum.
    s = np.nansum(np.sign(P_i) * np.sqrt(np.abs(P_i))) / (n * 1.0)

    return s


def spreadsheet_maker(df):
    """
    Compute several summary properties of EVERY column in this table, grouped by SOURCEID:
    - mean
    - median
    - std / rms
    - min
    - max
    - range (total, and 90/10)

    and also the following columns:
    - J_counts
    - H_counts
    - K_counts
    (possibly sorted by "good/noflag", "warn", "error" / from the ppErrBits)
    (see, for example, how Nicholas Cross's WSA databases handle this.)

    and some specifically variability-related columns:
    - chisq (J, H, K)
    - Stetson (JHK, JK, JH, HK)

    possibly with "robust" / "outlier-removed" versions of some or all of these.

    I'd like to perhaps group all these into a hierarchical data structure that looks like this:

    # ss : stands for SpreadSheet 
    spreadsheet.mean['JAPERMAG3']
    spread


    # JK! we're doing this entirely in pandas, no extra structs like a NamedTuple.
    spreadsheet['mean']['JAPERMAG3']

    """

    single_column_functions = [
        np.nanmean,
        np.nanmedian,
        np.nanmin,
        np.nanmax,
        np.nanstd,
        lambda x: np.nanmax(x) - np.nanmin(x),
        lambda x: np.nanpercentile(x, 90) - np.nanpercentile(x, 10),
    ]

    single_column_function_names = [
        "mean",
        "median",
        "min",
        "max",
        "std",
        "range",
        "range_9010",
    ]

    count_functions = [
        lambda x, y: np.sum(~np.isnan(x)),
        lambda x, y: np.sum(~np.isnan(x[y == 0])),
        lambda x, y: np.sum(~np.isnan(x[(0 < y) & (y < 256)])),
        lambda x, y: np.sum(~np.isnan(x[(256 <= y) & (y < 65536)])),
        lambda x, y: np.sum(~np.isnan(x[y >= 65536])),
    ]

    bands = ["J", "H", "K"]
    count_names = [
        lambda x: f"N_{x}",
        lambda x: f"N_{x}_good",
        lambda x: f"N_{x}_info",
        lambda x: f"N_{x}_warn",
        lambda x: f"N_{x}_severe",
    ]

    # rename this function.
    def f_mi(x):
        d = []
        primary_index = []
        secondary_index = []
        for fn, fn_name in zip(single_column_functions, single_column_function_names):
            for column in df.columns[1:]:

                d.append(fn(x[column]))
                primary_index.append(fn_name)
                secondary_index.append(column)

        for fn, fn_name in zip(count_functions, count_names):
            for band in bands:
                d.append(fn(x[f"{band}APERMAG3"], x[f"{band}PPERRBITS"]))
                primary_index.append("count")
                secondary_index.append(fn_name(band))

        return pd.Series(d, index=[primary_index, secondary_index])

    df_spreadsheet = df.groupby("SOURCEID").apply(f_mi)

    return df_spreadsheet


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
