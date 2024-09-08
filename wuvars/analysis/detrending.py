"""
Functions to remove secular trends in lightcurves in order to search for periodic variability.

"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from numpy.polynomial.polynomial import polyfit, polyval
from wuvars.analysis.periods import N_eval, f_max, f_min
from wuvars.plotting.lightcurve import (simple_lc_scatter_brokenaxes,
                                        simple_phased_lc_scatter_gridspec)


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

        try:
            fit_params = polyfit(times, mags, poly_order, w=1 / errs)
            fit_mag = polyval(date, fit_params)
        except TypeError:
            fit_mag = np.zeros_like(mags_raw)

        mag_mean = np.mean(mags)
        detrended_mags = mags_raw - fit_mag + mag_mean

        # overwrite

        detrended_dat[f"{band}APERMAG3"] = detrended_mags
        fit_dat[f"{band}APERMAG3"] = fit_mag

    return detrended_dat, fit_dat


def visualize_detrending(
    sid,
    dat,
    date_offset,
    start=0,
    end=np.inf,
    poly_order=1,
    nterms=1,
    cmap="jet_r",
    simple_breaks=None,
    plot_phase=False,
    force_period=False
):

    
    if len(cmap) == 2:
        cmap_1 = cmap[0]
        cmap_2 = cmap[1]
    else:
        cmap_1 = cmap
        cmap_2 = cmap

    fig_1 = simple_lc_scatter_brokenaxes(
        dat, sid, date_offset=date_offset, cmap=cmap_1, breaks=simple_breaks
    )

    dat_sid = dat.groups[dat.groups.keys["SOURCEID"] == sid]
    detrended, fit = poly_detrend(
        dat_sid, date_offset, start, end, poly_order=poly_order
    )

    fig_2 = simple_lc_scatter_brokenaxes(
        detrended.group_by("SOURCEID"), sid, cmap=cmap_2, breaks=[]
    )
    fig_3 = simple_lc_scatter_brokenaxes(
        fit.group_by("SOURCEID"), sid, cmap=cmap_2, breaks=[]
    )

    star_dat = detrended.group_by("SOURCEID")

    freq = np.linspace(f_min, f_max, N_eval)

    periods = []

    for band in ["J", "H", "K"]:
        # band = 'K'
        mask = star_dat[f"{band}APERMAG3"].mask
        times = star_dat["MEANMJDOBS"][~mask]
        mags = star_dat[f"{band}APERMAG3"][~mask]
        errs = star_dat[f"{band}APERMAG3ERR"][~mask]

        ls = LombScargle(times, mags, dy=errs, nterms=nterms)
        power = ls.power(freq, assume_regular_frequency=True).value

        min_freq = 1 / 100
        power[freq < min_freq] = 0
        power[np.abs(freq - 1) < 0.01] = 0

        # powermax = np.nanmax(power)
        fmax = freq[np.nanargmax(power)]
        #         fap = ls.false_alarm_probability(np.nanmax(power)).value

        plt.figure()
        plt.plot(1 / freq, power, "k", rasterized=True)
        plt.semilogx()
        plt.xlim(0.1, 10)

        print(f"{band} period: {1 / fmax:.4f}")
        try:
            amp = np.sqrt(np.sum(ls.model_parameters(fmax)[1:3] ** 2))
            print(f"{band} amp: {amp:.3f}")
        except:
            pass

        periods.append(1 / fmax)

    if plot_phase:

        if force_period:
            period = force_period
        else:

            if np.abs(periods[0] - periods[2])/periods[2] < 0.01:
                period = np.mean(periods)
                print("periods agree: ", periods)
            else:
                period = periods[2]
                print("periods disagree: ", periods, "; using K")

        fig_4 = simple_phased_lc_scatter_gridspec(
            star_dat, sid, period, cmap=cmap_2
        )
