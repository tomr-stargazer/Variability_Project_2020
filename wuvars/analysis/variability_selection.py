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

# let's write up the variability functions - I think I implemented these recently in
# stetson_2020.py. (in Variability_2019)


def delta(m, sigma_m):
    """
    Normalized residual / "relative error" for one observation.
    Used in Stetson's J variability index.

    INPUTS:
        m: a single magnitude measurement in a certain band
        sigma_m: the uncertainty on that measurement

    OUTPUTS:
        d: the "relative error"

    """

    n = m[~np.isnan(m)].size
    d = np.sqrt(n / (n - 1)) * (m - np.nanmean(m)) / sigma_m

    return d


# From c. 2012.
def S_threeband(j, sigma_j, h, sigma_h, k, sigma_k):
    """
    Computes the Stetson variability index J (renamed S) for one star that has
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

    n = np.sum(~np.isnan(j) & ~np.isnan(h) & ~np.isnan(k))

    # Perhaps hackish - but tests show it's needed to avoid np.inf when n==1.
    if n < 2:
        return np.nan

    d_j = np.sqrt(n / (n - 1)) * (j - np.nanmean(j)) / sigma_j
    d_h = np.sqrt(n / (n - 1)) * (h - np.nanmean(h)) / sigma_h
    d_k = np.sqrt(n / (n - 1)) * (k - np.nanmean(k)) / sigma_k

    P_i = np.array([d_j * d_h, d_h * d_k, d_j * d_k])

    # I originally had two sums going: one over P_i, and one over all the
    # elements of n, but then I realized that a single sum over all axes
    # did the exact same thing (I tested it) so now it's back to one sum.
    s = np.nansum(np.sign(P_i) * np.sqrt(np.abs(P_i))) / (n * 1.0)

    return s


def S_twoband(a, sigma_a, b, sigma_b):
    """
    Computes the Stetson variability index J for one star that has
    2 observations on each night. Uses Carpenter et al.'s notation.
    Simplified from the general expression assuming 2 observations every
      night.

    INPUTS:
        a: an array of magnitudes in arbitrary photometric band "A"
        sigma_a: an array of corresponding uncertainties
        b: an array of magnitudes in arbitrary photometric band "B"
        sigma_b: an array of corresponding uncertainties

    OUTPUTS:
        s: the Stetson variability index for 2 bands

    """

    n = np.sum(~np.isnan(a) & ~np.isnan(b))

    if n < 2:
        return np.nan

    d_a = np.sqrt(n / (n - 1)) * (a - np.nanmean(a)) / sigma_a
    d_b = np.sqrt(n / (n - 1)) * (b - np.nanmean(b)) / sigma_b

    P_i = d_a * d_b

    s = np.nansum(np.sign(P_i) * np.sqrt(np.abs(P_i))) / (n * 1.0)

    return s


def chisq(mag, err):
    """
    Difference between each value x_i and the mean of x, squared,
    divided by each error xerr_i, squared, summed for all i.

    https://en.wikipedia.org/wiki/Chi-squared_test#Pearson's_%CF%872_test

    """
    return np.nansum((mag - np.nanmean(mag)) ** 2 / err ** 2)


def reduced_chisq(mag, err):
    """
    Chisq divided by number of degrees of freedom (which is n-1).

    Caveat emptor: https://arxiv.org/abs/1012.3754

    """
    nu = np.sum(~np.isnan(mag)) - 1
    if nu == 0:
        return np.nan
    return 1 / nu * np.nansum((mag - np.nanmean(mag)) ** 2 / err ** 2)


def data_nuller(df, max_flags=256):
    """
    Applies rules to the dataframe to null out specific columns!

    Operates on the dataframe in-place.

    """

    null_val = -999999488.0

    # we care if:
    # - JAPERMAG is null, or
    # - JPPERRBITS > max_flags
    # THEN
    # - nan out JAPERMAG3, JAPERMAG3ERR, JMHPNT, JMHPNTERR, and JPPERRBITS

    j_nan = (df["JAPERMAG3"] == null_val) | (df["JPPERRBITS"] > max_flags)
    df.loc[j_nan, "JAPERMAG3"] = np.nan
    df.loc[j_nan, "JAPERMAG3ERR"] = np.nan
    df.loc[j_nan, "JMHPNT"] = np.nan
    df.loc[j_nan, "JMHPNTERR"] = np.nan
    df.loc[j_nan, "JPPERRBITS"] = np.nan

    # - HAPERMAG is null, or
    # - HPPERRBITS > max_flags
    # THEN
    # - nan out HAPERMAG3, HAPERMAG3ERR, JMHPNT, JMHPNTERR, HMKPNT, HMKPNTERR, and HPPERRBITS

    h_nan = (df["HAPERMAG3"] == null_val) | (df["HPPERRBITS"] > max_flags)
    df.loc[h_nan, "HAPERMAG3"] = np.nan
    df.loc[h_nan, "HAPERMAG3ERR"] = np.nan
    df.loc[h_nan, "JMHPNT"] = np.nan
    df.loc[h_nan, "JMHPNTERR"] = np.nan
    df.loc[h_nan, "HMKPNT"] = np.nan
    df.loc[h_nan, "HMKPNTERR"] = np.nan
    df.loc[h_nan, "HPPERRBITS"] = np.nan

    # - KAPERMAG is null
    # - KPPERRBITS > max_flags
    # THEN we wanna:
    # - nan out KAPERMAG3, KAPERMAG3ERR, HMKPNT, HMKPNTERR, and KPPERRBITS

    k_nan = (df["KAPERMAG3"] == null_val) | (df["KPPERRBITS"] > max_flags)
    df.loc[k_nan, "KAPERMAG3"] = np.nan
    df.loc[k_nan, "KAPERMAG3ERR"] = np.nan
    df.loc[k_nan, "HMKPNT"] = np.nan
    df.loc[k_nan, "HMKPNTERR"] = np.nan
    df.loc[k_nan, "KPPERRBITS"] = np.nan

    # raise an error if there's anything left by this
    assert ~np.any(df == null_val)


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

        # now we do the variability statistics!
        # one-band variability: just chisq on each band
        for band in bands:
            mag = x[f"{band}APERMAG3"]
            err = x[f"{band}APERMAG3ERR"]
            d.append(reduced_chisq(mag, err))
            primary_index.append("variability")
            secondary_index.append(band + "_red_chisq")

        # two-band variability: two-band Stetson
        for band in bands:
            # clever: iterate over permutation pairs of JHK by removing each
            band_pair = "JHK".replace(band, "")
            band_A, band_B = band_pair

            A_mag = x[f"{band_A}APERMAG3"]
            A_err = x[f"{band_A}APERMAG3ERR"]
            B_mag = x[f"{band_B}APERMAG3"]
            B_err = x[f"{band_B}APERMAG3ERR"]

            d.append(S_twoband(A_mag, A_err, B_mag, B_err))
            primary_index.append("variability")
            secondary_index.append(f"Stetson_{band_pair}")

        # three-band variability: three-band Stetson
        stetson_arguments = (
            x[f"{band}APERMAG3{err}"] for band in bands for err in ["", "ERR"]
        )
        d.append(S_threeband(*stetson_arguments))
        primary_index.append("variability")
        secondary_index.append("Stetson_JHK")

        return pd.Series(d, index=[primary_index, secondary_index])

    # # nullify
    # df[df == -999999488.0] = np.nan
    data_nuller(df)

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

    pass

    # v = None
    # ds = None
    # S_JHK = "Stetson_JHK"
    # S_JH = "Stetson_JH"
    # S_JK = "Stetson_JK"
    # S_HK = "Stetson_HK"

    # q2_all_indices = (
    #     (ds["count"]["N_J"] > 50)
    #     & (ds["count"]["N_J"] < 150)
    #     & (ds["count"]["N_H"] > 50)
    #     & (ds["count"]["N_H"] < 150)
    #     & (ds["count"]["N_K"] > 50)
    #     & (ds["count"]["N_K"] < 150)
    #     & (ds["max"]["JPPERRBITS"] == 0)
    #     & (ds["max"]["HPPERRBITS"] == 0)
    #     & (ds["max"]["KPPERRBITS"] == 0)
    #     & (ds["median"]["PSTAR"] > 0.75)
    # )

    # variable_indices = q2_all_indices & (
    #     (ds[v][S_JHK] > 2)
    #     | (ds[v][S_JH] > 2)
    #     | (ds[v][S_HK] > 2)
    #     | (ds[v][S_JK] > 2)
    #     | (ds[v]["J_red_chisq"] > 2)
    #     | (ds[v]["H_red_chisq"] > 2)
    #     | (ds[v]["K_red_chisq"] > 2)
    # )


# q1_old = (
#     (
#         (ds["count"]["N_J"] >= 50)
#         & (ds["count"]["N_J"] <= 125)
#         & (ds["mean"]["JAPERMAG3"] > 11)
#         & (ds["mean"]["JAPERMAG3"] < 17)
#         & (ds["count"]["N_J_info"] == 0)
#     )
#     | (
#         (ds["count"]["N_H"] >= 50)
#         & (ds["count"]["N_H"] <= 125)
#         & (ds["mean"]["HAPERMAG3"] > 11)
#         & (ds["mean"]["HAPERMAG3"] < 16)
#         & (ds["count"]["N_H_info"] == 0)
#     )
#     | (
#         (ds["count"]["N_K"] >= 50)
#         & (ds["count"]["N_K"] <= 125)
#         & (ds["mean"]["KAPERMAG3"] > 11)
#         & (ds["mean"]["KAPERMAG3"] < 16)
#         & (ds["count"]["N_K_info"] == 0)
#     )
# ) & (ds["median"]["PSTAR"] > 0.75)


# cand_case2 = (sp.pstar_median > 0.75) & (
#     (
#         (sp.N_j >= 50)
#         & (sp.N_j <= 125)
#         & (sp.j_mean > 11)  # J band criteria
#         & (sp.j_mean < 17)
#         & (sp.N_j_info == 0)
#     )
#     | (
#         (sp.N_h >= 50)
#         & (sp.N_h <= 125)
#         & (sp.h_mean > 11)  # H band criteria
#         & (sp.h_mean < 16)
#         & (sp.N_h_info == 0)
#     )
#     | (
#         (sp.N_k >= 50)
#         & (sp.N_k <= 125)
#         & (sp.k_mean > 11)  # K band criteria
#         & (sp.k_mean < 16)
#         & (sp.N_k_info == 0)
#     )
# )

# cand_case1 = (
#     (sp.pstar_median > 0.75)
#     & (
#         (sp.N_j >= 50)
#         & (sp.N_j <= 125)
#         & (sp.j_mean > 11)  # J band criteria
#         & (sp.j_mean < 17)
#         & (sp.N_j_info == 0)
#     )
#     & (
#         (sp.N_h >= 50)
#         & (sp.N_h <= 125)
#         & (sp.h_mean > 11)  # H band criteria
#         & (sp.h_mean < 16)
#         & (sp.N_h_info == 0)
#     )
#     & (
#         (sp.N_k >= 50)
#         & (sp.N_k <= 125)
#         & (sp.k_mean > 11)  # K band criteria
#         & (sp.k_mean < 16)
#         & (sp.N_k_info == 0)
#     )
# )
# q2_old = (
#     (
#         (ds["count"]["N_J"] >= 50)
#         & (ds["count"]["N_J"] <= 125)
#         & (ds["mean"]["JAPERMAG3"] > 11)
#         & (ds["mean"]["JAPERMAG3"] < 17)
#         & (ds["count"]["N_J_info"] == 0)
#     )
#     & (
#         (ds["count"]["N_H"] >= 50)
#         & (ds["count"]["N_H"] <= 125)
#         & (ds["mean"]["HAPERMAG3"] > 11)
#         & (ds["mean"]["HAPERMAG3"] < 16)
#         & (ds["count"]["N_H_info"] == 0)
#     )
#     & (
#         (ds["count"]["N_K"] >= 50)
#         & (ds["count"]["N_K"] <= 125)
#         & (ds["mean"]["KAPERMAG3"] > 11)
#         & (ds["mean"]["KAPERMAG3"] < 16)
#         & (ds["count"]["N_K_info"] == 0)
#     )
# ) & (ds["median"]["PSTAR"] > 0.75)


def q0_selection(ds):
    q0 = (
        (ds["count"]["N_J"] >= 50)
        | (ds["count"]["N_H"] >= 50)
        | (ds["count"]["N_K"] >= 50)
    )

    return q0


def q2_selection(ds):
    q2 = (
        (ds["count"]["N_J"] >= 50)
        & (ds["count"]["N_J"] < 150)
        & (ds["count"]["N_H"] >= 50)
        & (ds["count"]["N_H"] < 150)
        & (ds["count"]["N_K"] >= 50)
        & (ds["count"]["N_K"] < 150)
        & (ds["mean"]["JAPERMAG3"] > 11)
        & (ds["mean"]["HAPERMAG3"] > 11)
        & (ds["mean"]["KAPERMAG3"] > 11)
        & (ds["count"]["N_J"] == ds["count"]["N_J_good"])
        & (ds["count"]["N_H"] == ds["count"]["N_H_good"])
        & (ds["count"]["N_K"] == ds["count"]["N_K_good"])
        & (ds["median"]["PSTAR"] > 0.75)
    )
    return q2


def q1_selection(ds):
    q1 = (
        (
            (ds["count"]["N_J"] >= 50)
            & (ds["count"]["N_J"] < 125)
            & (ds["mean"]["JAPERMAG3"] > 11)
            & (ds["count"]["N_J"] == ds["count"]["N_J_good"])
        )
        | (
            (ds["count"]["N_H"] >= 50)
            & (ds["count"]["N_H"] < 125)
            & (ds["mean"]["HAPERMAG3"] > 11)
            & (ds["count"]["N_H"] == ds["count"]["N_H_good"])
        )
        | (
            (ds["count"]["N_K"] >= 50)
            & (ds["count"]["N_K"] < 125)
            & (ds["mean"]["KAPERMAG3"] > 11)
            & (ds["count"]["N_K"] == ds["count"]["N_K_good"])
        )
        & (ds["median"]["PSTAR"] > 0.75)
    )

    return q1


# autocandidates_old = sp.where(
#     (
#         (sp.N_j >= 50)
#         & (sp.N_j <= 125)
#         & (sp.j_mean > 11)  # J band criteria
#         & (sp.j_mean < 17)
#         & (sp.N_j_info == 0)  # J
#     )
#     | (  # J
#         (sp.N_h >= 50)
#         & (sp.N_h <= 125)
#         & (sp.h_mean > 11)  # H band criteria
#         & (sp.h_mean < 16)
#         & (sp.N_h_info == 0)  # H
#     )
#     | (  # H
#         (sp.N_k >= 50)
#         & (sp.N_k <= 125)
#         & (sp.k_mean > 11)  # K band criteria
#         & (sp.k_mean < 16)
#         & (sp.N_k_info == 0)  # K
#     )
# )  # K


def q1_variables(ds):
    q1v = (
        (
            (ds["count"]["N_J"] >= 50)
            & (ds["count"]["N_J"] < 150)
            & (ds["mean"]["JAPERMAG3"] > 11)
            & (ds["count"]["N_J"] == ds["count"]["N_J_good"])
            & (ds["variability"]["J_red_chisq"] > 2)
        )
        | (
            (ds["count"]["N_H"] >= 50)
            & (ds["count"]["N_H"] < 150)
            & (ds["mean"]["HAPERMAG3"] > 11)
            & (ds["count"]["N_H"] == ds["count"]["N_H_good"])
            & (ds["variability"]["H_red_chisq"] > 2)
        )
        | (
            (ds["count"]["N_K"] >= 50)
            & (ds["count"]["N_K"] < 150)
            & (ds["mean"]["KAPERMAG3"] > 11)
            & (ds["count"]["N_K"] == ds["count"]["N_K_good"])
            & (ds["variability"]["K_red_chisq"] > 2)
        )
        & (ds["median"]["PSTAR"] > 0.75)
    )

    return q1v
