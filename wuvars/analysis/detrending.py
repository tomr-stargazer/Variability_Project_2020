"""
Functions to remove secular trends in lightcurves in order to search for periodic variability.

"""

from numpy.polynomial.polynomial import polyfit, polyval
import numpy as np
import matplotlib.pyplot as plt


# from astropy.convolution import convolve, Box1DKernel
# rv_smooth = convolve(rv, Box1DKernel(100))
# rv_hpf = rv - rv_smooth


def poly_detrend(
    dat, date_offset=0, data_start=0, data_end=np.inf, poly_order=1, **kwargs
):
    """
    Returns a dataset for a single star that has been detrended.

    """
    # dat has already been through a step like this:
    # dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]
    # dg: astropy table that has been grouped by SOURCEID

    unpruned_date = dat["MEANMJDOBS"] - date_offset

    season_selection = (unpruned_date >= data_start) & (unpruned_date <= data_end)

    dat_pruned = dat[season_selection]
    dp = dat_pruned

    detrended_dat = dp.copy()
    fit_dat = dp.copy()

    # generate the polyfit
    bands = ["J", "H", "K"]

    date = dp["MEANMJDOBS"] - date_offset

    for band in bands:

        mask = dp[f"{band}APERMAG3"].mask
        times = dp["MEANMJDOBS"][~mask] - date_offset

        mags_raw = dp[f"{band}APERMAG3"]
        errs_raw = dp[f"{band}APERMAG3ERR"]

        mags = dp[f"{band}APERMAG3"][~mask]
        errs = dp[f"{band}APERMAG3ERR"][~mask]

        fit_params = polyfit(times, mags, poly_order, w=1 / errs)
        fit_mag = polyval(date, fit_params)

        mag_mean = np.mean(mags)
        detrended_mags = mags_raw - fit_mag + mag_mean

        # overwrite

        detrended_dat[f"{band}APERMAG3"] = detrended_mags
        fit_dat[f"{band}APERMAG3"] = fit_mag

    return detrended_dat, fit_dat
