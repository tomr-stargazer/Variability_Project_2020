"""
I'm going to calculate Q, the quasiperiodicity metric identified by Cody et al. 2014.

I'm adopting the simplification described in Bredall+20 and Hillenbrand+22.

"""

import numpy as np
from astropy.timeseries import LombScargle
from wuvars.analysis.periods import N_eval, f_max, f_min, freq

# from: Bredall et al. 2020
# at: https://academic.oup.com/mnras/article/496/3/3257/5855495?login=false
# "To produce the residual curve, we create a window function (Dirac comb) out of
# the original light curve, i.e. values of 1 when an observation takes place and values
# of 0 elsewhere. We utilize astropy lombscargle (VanderPlas et al. 2012;
# VanderPlas & Ivezić 2015) to find the periodograms of both the original curve and
# the Dirac comb window function. We find all frequencies in the window function with a
# power greater than 2 standard deviations above the mean and ignore any peaks in the
# periodogram of the light curve within 0.03 Hz of the peaks in the window function.
# The curve is then phased using the remaining period with the highest power. We
# boxcar-smooth this folded light curve with a window size of 25 per cent of the period
# and subtract the smoothed model from the folded curve. We define σresid as the
# standard deviation of this residual curve."

# FOR EACH band:
# do a period run on the signal
# do a period run on the dirac comb
# pick the strongest period that isn't excluded by the windowing function
#   a. compute the window function

#     from https://github.com/jakevdp/PracticalLombScargle/blob/master/figures/LINEAR_binary.ipynb
#     """from astropy.timeseries import LombScargle
#     ls_window = LombScargle(data.t, 1, fit_mean=False, center_data=False)
#     freq, power = ls_window.autopower(minimum_frequency=0.001,
#                                       maximum_frequency=10)

#   b. compute the periodogram for the data

#   c. select out the

# fig, ax = plt.subplots(figsize=(8, 3))
# ax.semilogx(1. / freq, power, '-k');"""


def smooth_phase(t, mag, period):
    """
    Given a period and a time series signal, I want the "smoothed" phase-folded curve.

    """
    # window_width
    ww = 0.25

    # setup
    phased_t = (t % period) / period
    smoothed_mag = np.copy(mag)

    for i in range(len(phased_t)):
        ti = phased_t[i]
        in_window = (np.abs(ti - phased_t) < ww / 2) | (
            np.abs(ti - phased_t) > (1 - ww / 2)
        )
        smoothed_mag[i] = np.mean(mag[in_window])

    return smoothed_mag
    # the residual should actually be literally the original mag minus the smoothed_mag


def Q_score(mags, err, residual):
    """ Calculates Q given mags, errs, and residuals. """

    rms_resid = np.std(residual)
    rms_raw = np.std(mags)

    Q = (rms_resid ** 2 - err ** 2) / (rms_raw ** 2 - err ** 2)

    return Q


# boxcar smooth it (rolling average) with a window = 0.25 of the period
# create a residual lightcurve by subtracting (with interpolation) the smooth one from the real one
# calculate Q via the equation


def compute_Q_automatically(dataset, sid, **kwargs):
    """ 
    Given a photometry dataset and a SID, I want you to compute Q.

    This function does the period selection (for the residual) automatically.

    THis may be slightly inefficient (probably it would be better to not pass the whole dataset
    into the function on every iterative step) but I'm just prototyping here.

    The more "mature" version will probably be a script that loads and references
    the photometric dataset just once, and then loops through all the SIDs to calculate Q.
    But for now this is fine.
    """

    dg = dataset
    dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

    Q_dict = {}

    # FOR EACH band:
    for band in ["J", "H", "K"]:
        mask = dat[f"{band}APERMAG3"].mask
        times = dat["MEANMJDOBS"][~mask]
        mags = dat[f"{band}APERMAG3"][~mask]
        errs = dat[f"{band}APERMAG3ERR"][~mask]

        try:
            Q_dict[band] = compute_Q_band(times, mags, errs, **kwargs)
        except ValueError:
            Q_dict[band] = (np.nan, np.nan)

    return Q_dict


def compute_Q_band(times, mags, errs, verbose=True):

    # do a period run on the dirac comb
    ls_window = LombScargle(times, 1, fit_mean=False, center_data=False)
    power_window = ls_window.power(freq)

    # do a period run on the signal
    ls = LombScargle(times, mags, dy=errs)
    power = ls.power(freq, assume_regular_frequency=True)

    # pick the strongest period that isn't excluded by the windowing function
    exclude = power_window > 2 * np.nanstd(power_window)

    fmax = freq[~exclude][np.nanargmax(power[~exclude])]
    if verbose:
        print(f"1/fmax = {1/fmax:.3f}")
    smooth_mag = smooth_phase(times, mags, 1 / fmax)
    residual = mags - smooth_mag

    sigma = np.median(errs)

    Q = Q_score(mags, sigma, residual)

    return Q, fmax
