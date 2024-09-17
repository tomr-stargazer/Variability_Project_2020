"""
I wanted to pull out the period results by themselves.


"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
import wuvars.analysis.variability_selection as sv
from numpy.polynomial.polynomial import polyfit, polyval
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.load_periodics_v4 import (ic_periods, ngc_periods,
                                               select_periodic_variables_v4)
from wuvars.analysis.spectral_type_to_number import get_SpT_from_num
from wuvars.analysis.variability_selection_curved import (curve_Stetson, sv_hk,
                                                          sv_jh, sv_jhk, sv_jk)
from wuvars.data import photometry, quality_classes, spreadsheet

warnings.filterwarnings("ignore")

# This part pulls up the source information from our position-matching scheme.
ngc_match = match_ngc()
ic_match = match_ic()

names = ["ngc", "ic"]
wserv_dict = {"ngc": 7, "ic": 8}
fullname_dict = {"ngc": "NGC 1333", "ic": "IC 348"}
match_dict = {"ngc": ngc_match, "ic": ic_match}
spread_dict = {"ngc": spreadsheet.load_wserv_v2(7), "ic": spreadsheet.load_wserv_v2(8)}
q_dict = {"ngc": quality_classes.load_q(7), "ic": quality_classes.load_q(8)}
per_dict = {"ngc": ngc_periods, "ic": ic_periods}

for name in names:

    print("=========")
    print(f"{fullname_dict[name]}:")

    periods = per_dict[name]
    match = match_dict[name]
    spread = spread_dict[name]

    # let's divide up by infrared excess
    ir_exc = match.approved["IRexc"] == "yes"

    fig_peramp = plt.figure(figsize=(10, 10))
    ax_dummy = fig_peramp.add_subplot(111)
    ax_dummy.axis("off")

    ax_peramp = fig_peramp.add_subplot(223)
    ax_per_hist = fig_peramp.add_subplot(221, sharex=ax_peramp)
    ax_amp_hist = fig_peramp.add_subplot(224, sharey=ax_peramp)

    # So. What do we want on this figure?
    # Because... I kind of want...

    per_nir_50 = np.nanmedian(periods["Period"][~ir_exc])
    per_nir_84 = np.nanpercentile(periods["Period"][~ir_exc], 84)
    per_nir_16 = np.nanpercentile(periods["Period"][~ir_exc], 16)

    per_ir_50 = np.nanmedian(periods["Period"][ir_exc])
    per_ir_84 = np.nanpercentile(periods["Period"][ir_exc], 84)
    per_ir_16 = np.nanpercentile(periods["Period"][ir_exc], 16)

    amp_nir_50 = np.nanmedian(periods["Amp"][~ir_exc])
    amp_nir_84 = np.nanpercentile(periods["Amp"][~ir_exc], 84)
    amp_nir_16 = np.nanpercentile(periods["Amp"][~ir_exc], 16)

    amp_ir_50 = np.nanmedian(periods["Amp"][ir_exc])
    amp_ir_84 = np.nanpercentile(periods["Amp"][ir_exc], 84)
    amp_ir_16 = np.nanpercentile(periods["Amp"][ir_exc], 16)

    print(
        f"Per nir 50 +/- 1 sigma: {per_nir_50:.1f} +{per_nir_84- per_nir_50:.1f} -{per_nir_50- per_nir_16:.1f}"
    )

    print(
        f"Per IR 50 +/- 1 sigma: {per_ir_50:.1f} +{per_ir_84- per_ir_50:.1f} -{per_ir_50- per_ir_16:.1f}"
    )

    print("--")
    print(
        f"Amp nir 50 +/- 1 sigma: {amp_nir_50:.1e} +{amp_nir_84- amp_nir_50:.1e} -{amp_nir_50- amp_nir_16:.1e}"
    )

    print(
        f"Amp IR 50 +/- 1 sigma: {amp_ir_50:.1e} +{amp_ir_84- amp_ir_50:.1e} -{amp_ir_50- amp_ir_16:.1e}"
    )

    ax_peramp.plot(
        periods["Period"][ir_exc], periods["Amp"][ir_exc], "r.", label="IR excess"
    )
    ax_peramp.plot(
        periods["Period"][~ir_exc], periods["Amp"][~ir_exc], "k.", label="no IR excess"
    )
    ax_peramp.loglog()
    ax_peramp.set_xlabel("Period (d)")
    ax_peramp.set_ylabel("Amplitude (mag)")
    ax_peramp.legend()

    min_per_bin = np.log10(np.nanmin(periods["Period"]))
    max_per_bin = np.log10(np.nanmax(periods["Period"]))
    print("Log Per bins: ", min_per_bin, max_per_bin)

    min_amp_bin = np.log10(np.nanmin(periods["Amp"]))
    max_amp_bin = np.log10(np.nanmax(periods["Amp"]))
    print("Log amp bins: ", min_amp_bin, max_amp_bin)

    min_per_bin = min(-0.11314264643119235, -0.2091950052836405) * 0.99
    max_per_bin = max(1.744725685332582, 1.2013475449906188) * 1.01

    min_amp_bin = min(-2.3159896476462887, -2.091253337507326) * 0.99
    max_amp_bin = max(-0.6408834417764375, -0.9407821623586071) * 1.01

    ax_per_hist.hist(
        periods["Period"][ir_exc],
        bins=np.logspace(min_per_bin, max_per_bin, 10),
        edgecolor="r",
        facecolor="None",
        hatch="//",
        histtype="stepfilled",
        label="IR excess",
    )
    ax_per_hist.hist(
        periods["Period"][~ir_exc],
        bins=np.logspace(min_per_bin, max_per_bin, 10),
        edgecolor="k",
        facecolor="None",
        hatch="..",
        histtype="stepfilled",
        label="no IR excess",
    )

    ax_per_hist.axvline(per_ir_50, color="r", linestyle="--")
    ax_per_hist.axvline(per_nir_50, color="k", linestyle=":")

    ax_amp_hist.hist(
        periods["Amp"][ir_exc],
        bins=np.logspace(min_amp_bin, max_amp_bin, 10),
        orientation="horizontal",
        edgecolor="r",
        facecolor="None",
        hatch="//",
        histtype="stepfilled",
        label="IR excess",
    )

    ax_amp_hist.hist(
        periods["Amp"][~ir_exc],
        bins=np.logspace(min_amp_bin, max_amp_bin, 10),
        orientation="horizontal",
        edgecolor="k",
        facecolor="None",
        hatch="..",
        histtype="stepfilled",
        label="no IR excess",
    )

    ax_amp_hist.axhline(amp_ir_50, color="r", linestyle="--")
    ax_amp_hist.axhline(amp_nir_50, color="k", linestyle=":")

    ax_amp_hist.legend()
    ax_per_hist.legend()

    ax_peramp.axhline(amp_ir_50, color="r", linestyle="--", xmax=2, clip_on=False)
    ax_peramp.axhline(amp_nir_50, color="k", linestyle=":", xmax=2, clip_on=False)
    ax_peramp.axvline(per_ir_50, color="r", linestyle="--", ymax=2, clip_on=False)
    ax_peramp.axvline(per_nir_50, color="k", linestyle=":", ymax=2, clip_on=False)

    ax_dummy.set_title(
        f"Period-Amplitude distribution for periodic variables in {fullname_dict[name]},\n"
        "disaggregated by whether they have an IR excess"
    )

    # for each star forming region
    # x: period
    # y: amplitude
    # black dots: no IR Exc
    # red dots: IR exc
    # x axis histogram for each
    # y axis histogram for each

    # THEN

    # fit a line
    xs_ir = periods["SpT"][ir_exc]
    ys_ir = periods["Amp"][ir_exc]
    good_ir = np.isfinite(ys_ir)
    plot_xs = np.arange(0, 13)

    fit_params_ir = polyfit(xs_ir[good_ir], ys_ir[good_ir], 1)
    fit_ys_ir = polyval(plot_xs, fit_params_ir)

    xs_nir = periods["SpT"][~ir_exc]
    ys_nir = periods["Amp"][~ir_exc]
    good_nir = np.isfinite(ys_nir)

    fit_params_nir = polyfit(xs_nir[good_nir], ys_nir[good_nir], 1)
    fit_ys_nir = polyval(plot_xs, fit_params_nir)

    fig_sptamp = plt.figure(figsize=(6, 6))
    ax_sptamp = fig_sptamp.add_subplot(111)

    ax_sptamp.plot(
        periods["SpT"][ir_exc], periods["Amp"][ir_exc], "r.", label="IR excess"
    )
    ax_sptamp.plot(plot_xs, fit_ys_ir, "r--")

    ax_sptamp.plot(
        periods["SpT"][~ir_exc], periods["Amp"][~ir_exc], "k.", label="no IR excess"
    )
    ax_sptamp.plot(plot_xs, fit_ys_nir, "k:")

    ax_sptamp.semilogy()
    ax_sptamp.set_xlabel("SpT")
    ax_sptamp.set_ylabel("Amplitude (mag)")
    ax_sptamp.legend()

    # x: spectral type
    # y: amplitude
    # (rest is same as above)
    # Check to see if (a) there's a bulk offset but primarily (b) is there a trend with
    # amplitude versus spectral type for either population?

    ys_per_ir = periods["Period"][ir_exc]
    good_per_ir = np.isfinite(ys_per_ir)
    plot_xs = np.arange(0, 13)

    fit_params_per_ir = polyfit(xs_ir[good_per_ir], ys_per_ir[good_per_ir], 1)
    fit_ys_per_ir = polyval(plot_xs, fit_params_per_ir)

    ys_per_nir = periods["Period"][~ir_exc]
    good_per_nir = np.isfinite(ys_per_nir)

    fit_params_per_nir = polyfit(xs_nir[good_per_nir], ys_per_nir[good_per_nir], 1)
    fit_ys_per_nir = polyval(plot_xs, fit_params_per_nir)

    fig_sptper = plt.figure(figsize=(6, 6))

    ax_sptper = fig_sptper.add_subplot(111)

    ax_sptper.plot(
        periods["SpT"][ir_exc], periods["Period"][ir_exc], "r.", label="IR excess"
    )
    ax_sptper.plot(plot_xs, fit_ys_per_ir, "r--")

    ax_sptper.plot(
        periods["SpT"][~ir_exc], periods["Period"][~ir_exc], "k.", label="no IR excess"
    )
    ax_sptper.plot(plot_xs, fit_ys_per_nir, "k:")

    ax_sptper.semilogy()
    ax_sptper.set_xlabel("SpT")
    ax_sptper.set_ylabel("Period (d)")
    ax_sptper.legend()
    # x: spectral type
    # y: period
    # (rest is same as above)
    # Check to see if (a) there's a bulk offset but primarily (b) is there a trend with
    # period versus spectral type for either population?
    # IN PARTICULAR! periodic non IR exc sources should be strongly regarded as rotators
    # and therefore any trend with non-IR exc periods vs spectral type is a
    # mass-rotation sequence

    fig_pers_binned = plt.figure(figsize=(6, 6))

    ax_pers_binned = fig_pers_binned.add_subplot(111)

    bin_edges = [0, 2, 4, 6, 8, 15]
    left_edges = bin_edges[:-1]
    right_edges = bin_edges[1:]

    for left, right in zip(left_edges, right_edges):

        for ir_status, color, offset in zip([ir_exc, ~ir_exc], ["r", "k"], [0.2, 0]):

            in_this_bin = (match.approved["SpT"] >= left) & (
                match.approved["SpT"] < right
            )
            print(left, right)

            # collect periods

            median_period = np.nanmedian(periods["Period"][ir_status & in_this_bin])
            per_84 = np.nanpercentile(periods["Period"][ir_status & in_this_bin], 84)
            per_16 = np.nanpercentile(periods["Period"][ir_status & in_this_bin], 16)

            up_err = per_84 - median_period
            down_err = median_period - per_16

            ax_pers_binned.errorbar(
                [(left + right) / 2 + offset],
                [median_period],
                fmt=".",
                color=color,
                yerr=np.array([[down_err, up_err]]).T,
            )
            ax_pers_binned.hlines(y=median_period, xmin=left, xmax=right, color=color)

            center_dict = {
                "verticalalignment": "center",
                "horizontalalignment": "center",
                "color": color,
            }
            ax_pers_binned.text(left, median_period, "[", center_dict)
            ax_pers_binned.text(right, median_period, ")", center_dict)

    ax_pers_binned.semilogy()
    ax_pers_binned.set_xlabel("SpT")
    ax_pers_binned.set_ylabel("Period (d)")

    ###

    fig_amps_binned = plt.figure(figsize=(6, 6))

    ax_amps_binned = fig_amps_binned.add_subplot(111)

    bin_edges = [0, 2, 4, 6, 8, 15]
    left_edges = bin_edges[:-1]
    right_edges = bin_edges[1:]

    for left, right in zip(left_edges, right_edges):

        for ir_status, color, offset in zip([ir_exc, ~ir_exc], ["r", "k"], [0.2, 0]):

            in_this_bin = (match.approved["SpT"] >= left) & (
                match.approved["SpT"] < right
            )
            print(left, right)

            # collect periods

            median_amp = np.nanmedian(periods["Amp"][ir_status & in_this_bin])
            amp_84 = np.nanpercentile(periods["Amp"][ir_status & in_this_bin], 84)
            amp_16 = np.nanpercentile(periods["Amp"][ir_status & in_this_bin], 16)

            up_err = amp_84 - median_amp
            down_err = median_amp - amp_16

            ax_amps_binned.errorbar(
                [(left + right) / 2 + offset],
                [median_amp],
                fmt=".",
                color=color,
                yerr=np.array([[down_err, up_err]]).T,
            )
            ax_amps_binned.hlines(y=median_amp, xmin=left, xmax=right, color=color)

            center_dict = {
                "verticalalignment": "center",
                "horizontalalignment": "center",
                "color": color,
            }
            ax_amps_binned.text(left, median_amp, "[", center_dict)
            ax_amps_binned.text(right, median_amp, ")", center_dict)

    ax_amps_binned.semilogy()
    ax_amps_binned.set_xlabel("SpT")
    ax_amps_binned.set_ylabel("Amplitude (mag)")        

    # okay. I want to do some
