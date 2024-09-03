"""
This is a script that will attempt to remove the "final" outliers in 
graded-clipped-scrubbed-whatever data.

Here's what I wrote in `immediate post cool stars notes`:

> I’m going to resolve the issue in NGC 1333 of conspicuous “outlier” 
photometric datapoints. 
> I have a hunch that some combination of “grade” filtering plus 
“check if this data point is more than three times the photometric uncertainty 
from its neighbors plus if (e.g.) H and K are also not deviant” will handle it. 
I believe this is scientifically justifiable even if it is slightly ad hoc - 
it is, essentially, me formalizing the rule I already use to “by eye” judge whether 
there is any reliable physical meaning in any given “large outlier” data point. 
It is not simply “removing data that I don’t like”.
> I’ll also check to see if IC 348 needs a similar treatment.


"""

import pdb

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from brokenaxes import brokenaxes
from wuvars.plotting.lightcurve_helpers import orion_cmap, produce_xlims


# To some extent this belongs in wuvars.plotting but - whatever, it's here now.
def simple_lc_scatter_brokenaxes_grade(
    dg, sid, date_offset=None, pad=5, xlims=None, breaks=None, **kwargs
):
    # dg: astropy table that has been grouped by SOURCEID

    dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

    # set up data

    if date_offset is None:
        date_offset = np.floor(np.min(dat["MEANMJDOBS"]))
    date = dat["MEANMJDOBS"] - date_offset

    j = dat["JAPERMAG3"]
    h = dat["HAPERMAG3"]
    k = dat["KAPERMAG3"]

    j_e = dat["JAPERMAG3ERR"]
    h_e = dat["HAPERMAG3ERR"]
    k_e = dat["KAPERMAG3ERR"]

    j_grade = dat["JGRADE"]
    h_grade = dat["HGRADE"]
    k_grade = dat["KGRADE"]

    # set up plot
    # OUR brokenaxes changes will all go here
    # need:
    #   - how much MJD to subtract
    #   - where to put the breaks (+how many breaks)
    #   - ... that's p much it, right?

    fig = plt.figure(figsize=(10, 6), dpi=80, facecolor="w", edgecolor="k")

    gs0 = fig.add_gridspec(1, 2, width_ratios=[6, 3], wspace=0.2)
    gs_left = gs0[0].subgridspec(3, 1, hspace=0)
    gs_right = gs0[1].subgridspec(1, 2, width_ratios=(2, 0.075), wspace=0.15)
    gs01 = gs_right[0].subgridspec(2, 1)

    bax_kwargs = dict(despine=False, d=0.0075, tilt=60, wspace=0.04)

    # xlims = (
    #     (56848 - date_offset, 56900 - date_offset),
    #     (56920 - date_offset, 57080 - date_offset),
    #     (57200 - date_offset, 57250 - date_offset),
    # )

    # xlims = ((0, 30), (80, 225), (370,390))
    # xlims = [(0.0, 21.0), (84.0, 222.0), (377.0, 385.0)]
    if xlims is None:
        if breaks is None:
            # baked-in default since I was prototyping on IC348; possibly a bad default
            breaks = [50, 350]
        xlims = produce_xlims(date, breaks=breaks, pad=pad)
        # print(xlims)

    ax_j = brokenaxes(xlims=xlims, subplot_spec=gs_left[0, 0], **bax_kwargs)
    ax_h = brokenaxes(xlims=xlims, subplot_spec=gs_left[1, 0], **bax_kwargs)
    ax_k = brokenaxes(xlims=xlims, subplot_spec=gs_left[2, 0], **bax_kwargs)

    # Make the broken zones gray. (Uses some details from  )
    ax_j.big_ax.set_zorder(-100)
    ax_h.big_ax.set_zorder(-100)
    ax_k.big_ax.set_zorder(-100)
    ax_j.big_ax.set_facecolor("0.9")
    ax_h.big_ax.set_facecolor("0.9")
    ax_k.big_ax.set_facecolor("0.9")

    ax_khk = fig.add_subplot(gs01[0, 0])
    ax_jhk = fig.add_subplot(gs01[1, 0])

    # ax_cbar = fig.add_subplot(gs_right[1])

    ax_j.tick_params(labelbottom=False)
    ax_h.tick_params(labelbottom=False)

    # bottom = 0.1
    # height = 0.25
    # left = 0.075
    # width = 0.5

    # ax_k = fig.add_axes((left, bottom, width, height))
    # ax_h = fig.add_axes((left, bottom + 0.3, width, height), sharex=ax_k)
    # ax_j = fig.add_axes((left, bottom + 0.6, width, height), sharex=ax_k)

    # ax_jhk = fig.add_axes((0.65, bottom, 0.3, 0.375))
    # ax_khk = fig.add_axes((0.65, bottom + 0.475, 0.3, 0.375))

    fig.ax_k = ax_k
    fig.ax_j = ax_j
    fig.ax_h = ax_h
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    sc_j = ax_j.scatter(
        date,
        j,
        c=j_grade,
        vmin=None,
        vmax=None,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )
    sc_h = ax_h.scatter(
        date,
        h,
        c=h_grade,
        vmin=None,
        vmax=None,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )
    sc_k = ax_k.scatter(
        date,
        k,
        c=k_grade,
        vmin=None,
        vmax=None,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )

    ax_j.errorbar(
        date,
        j,
        yerr=j_e,
        fmt="None",
        ecolor="k",
        ms=2,
        elinewidth=0.5,
        zorder=-1,
        alpha=0.5,
    )
    ax_h.errorbar(
        date,
        h,
        yerr=h_e,
        fmt="None",
        ecolor="k",
        ms=2,
        elinewidth=0.5,
        zorder=-1,
        alpha=0.5,
    )
    ax_k.errorbar(
        date,
        k,
        yerr=k_e,
        fmt="None",
        ecolor="k",
        ms=2,
        elinewidth=0.5,
        zorder=-1,
        alpha=0.5,
    )

    sc_jhk = ax_jhk.scatter(
        h - k,
        j - h,
        c=np.min((j_grade, h_grade, k_grade), axis=0),
        vmin=None,
        vmax=None,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )

    sc_khk = ax_khk.scatter(
        h - k,
        k,
        c=np.min((h_grade, k_grade), axis=0),
        vmin=None,
        vmax=None,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )

    ax_jhk.errorbar(
        h - k,
        j - h,
        xerr=(h_e ** 2 + k_e ** 2) ** 0.5,
        yerr=(h_e ** 2 + j_e ** 2) ** 0.5,
        fmt="None",
        ecolor="k",
        ms=2,
        elinewidth=0.5,
        zorder=-1,
        alpha=0.1,
    )

    ax_khk.errorbar(
        h - k,
        k,
        xerr=(h_e ** 2 + k_e ** 2) ** 0.5,
        yerr=k_e,
        fmt="None",
        ecolor="k",
        ms=2,
        elinewidth=0.5,
        zorder=-1,
        alpha=0.1,
    )

    # cbar = fig.colorbar(
    #     sc_j[0], cax=ax_cbar
    # )  # This should really be changed to the method

    ax_j.invert_yaxis()
    ax_h.invert_yaxis()
    ax_k.invert_yaxis()
    ax_khk.invert_yaxis()

    ax_j.set_ylabel("J", labelpad=40, fontdict={"rotation": "horizontal"})
    ax_h.set_ylabel("H", labelpad=40, fontdict={"rotation": "horizontal"})
    ax_k.set_ylabel("K", labelpad=40, fontdict={"rotation": "horizontal"})

    # \u2212 is a proper minus sign (better than the hyphen character `-`)
    ax_k.set_xlabel(f"MJD \u2212 {date_offset}", labelpad=20)
    # cbar.set_label(f"MJD \u2212 {date_offset}")

    ax_jhk.set_xlabel("H-K")
    ax_jhk.set_ylabel("J-H")  # , {'rotation':'horizontal'})
    ax_khk.set_xlabel("H-K")
    ax_khk.set_ylabel("K")  # , {'rotation':'horizontal'})

    plt.colorbar(sc_jhk)

    return fig


def calculate_diff(dat, sid, band="J", verbose=False):

    dt = dat.groups[dat.groups.keys["SOURCEID"] == sid]

    # dt

    j_vals = dt["JAPERMAG3"]
    j_sigs = dt["JAPERMAG3ERR"]

    diff = np.zeros_like(j_vals)

    l = len(j_vals)

    for i in range(l):

        if verbose:
            print("i=", i)

        left_diff = np.nan
        right_diff = np.nan
        diff[i] = np.nan

        if not np.ma.is_masked(j_vals[i]):

            a_l = 1

            while (i - a_l) >= 0:

                left_diff = (j_vals[i - a_l] - j_vals[i]) / j_sigs[i]

                if verbose:
                    print("Left diff:", left_diff, "a_l = ", a_l)

                if i == 8:
                    pdb.set_trace()

                if not np.ma.is_masked(left_diff):
                    if verbose:
                        print("  Breaking -  we accept left_diff")
                    break
                else:
                    if verbose:
                        print("  Not breaking - we don't accept left_diff")
                    a_l += 1
                    if verbose:
                        print("  a_l = ", a_l)
            if verbose:
                print("a_l: ", a_l, "i-a_l=", i - a_l)
                print("")

            a_r = 1

            while (i + a_r) < l:

                right_diff = (j_vals[i + a_r] - j_vals[i]) / j_sigs[i]

                if not np.ma.is_masked(right_diff):
                    print("Breaking")
                    break
                else:
                    print("Not breaking")
                    a_r += 1
            if verbose:
                print("a_r: ", a_r, "i+a_r", i + a_r, "l=", l)

                print(left_diff, right_diff)

            try:
                diffs = np.array([left_diff, right_diff])
                abs_diffs = np.array([np.abs(left_diff), np.abs(right_diff)])
                diff_abs = min(abs_diffs)
                diff_sgn = np.sign(diffs[diff_abs == abs_diffs])[0]

                diff[i] = diff_sgn * diff_abs
            except IndexError:
                diff[i] = np.nan
        else:
            if verbose:
                print(f"i={i}: J={j_vals[i]} is masked")

    #         if i>3:
    #             break
    return diff


def calculate_diffs(dat, sid):
    """
    This calculates a `diff` metric which essentially describes
    how far each data point is from its neighbors, scaled by the photometric error.

    Most of the logic in the code deals with situations in which the nearest neighbor
    data point is `nan` or otherwise invalid, prompting the code to progressively 
    search further right (or left) until it finds valid data to compare against.
    """

    dt = dat.groups[dat.groups.keys["SOURCEID"] == sid]

    diffs_list = []

    for band in ["J", "H", "K"]:

        b_vals = dt[f"{band}APERMAG3"]
        b_sigs = dt[f"{band}APERMAG3ERR"]

        diff = np.zeros_like(b_vals)

        len_ = len(b_vals)

        for i in range(len_):

            left_diff = np.nan
            right_diff = np.nan
            diff[i] = np.nan

            if not np.ma.is_masked(b_vals[i]):

                a_l = 1

                while (i - a_l) >= 0:

                    left_diff = (b_vals[i - a_l] - b_vals[i]) / b_sigs[i]

                    if not np.ma.is_masked(left_diff):
                        break
                    else:
                        a_l += 1

                a_r = 1

                while (i + a_r) < len_:

                    right_diff = (b_vals[i + a_r] - b_vals[i]) / b_sigs[i]

                    if not np.ma.is_masked(right_diff):
                        break
                    else:
                        a_r += 1

                try:
                    diffs = np.array([left_diff, right_diff])
                    abs_diffs = np.array([np.abs(left_diff), np.abs(right_diff)])
                    diff_abs = min(abs_diffs)
                    diff_sgn = np.sign(diffs[diff_abs == abs_diffs])[0]

                    diff[i] = diff_sgn * diff_abs
                except IndexError:
                    diff[i] = np.nan

        diffs_list.append(diff)
    return diffs_list


def vprint(*args, verbose=True, **kwargs):
    if verbose:
        print(*args, **kwargs)


def identify_outliers(dat, sid, diffs, date_offset=56141, verbose=False):

    j_diff, h_diff, k_diff = diffs

    siddat_ = dat.groups[dat.groups.keys["SOURCEID"] == sid]

    vprint(f"For source {sid}: ", verbose=verbose)

    j_diff_mad = scipy.stats.median_abs_deviation(j_diff, nan_policy="omit")
    h_diff_mad = scipy.stats.median_abs_deviation(h_diff, nan_policy="omit")
    k_diff_mad = scipy.stats.median_abs_deviation(k_diff, nan_policy="omit")

    k_mad = 1.4826
    bands = ["J", "H", "K"]
    diffs = {"J": j_diff, "H": h_diff, "K": k_diff}
    mads = {"J": j_diff_mad, "H": h_diff_mad, "K": k_diff_mad}

    results = {"J": [], "H": [], "K": []}

    for i in range(len(j_diff)):

        for b in range(len(bands)):
            b0 = bands[b]  # J, H, K
            b1 = bands[b - 1]  # K, J, H
            b2 = bands[b - 2]  # H, K, J

            this_datum = diffs[b0][i]
            this_grade = siddat_[f"{b0}GRADE"][i]
            this_date = siddat_["MEANMJDOBS"][i]

            if (np.abs(this_datum) > max(3 * k_mad * mads[b0], 5)) and (
                this_grade < 0.98
            ):
                vprint(
                    f"Found an outlier! {b0}: {this_datum} and {this_grade}",
                    verbose=verbose,
                )
                vprint(
                    f"  at date={this_date-date_offset:.1f} (i={i})", verbose=verbose
                )
                #             vprint(f"  Found an outlier! {this_datum=} and {this_grade=}")

                # pdb.set_trace()

                if (np.abs(diffs[b1][i]) > 3 * k_mad * mads[b1]) or (
                    np.abs(diffs[b2][i]) > 3 * k_mad * mads[b2]
                ):
                    vprint(
                        "*** But there are outliers in other bands at that time too. No flag. ***",
                        verbose=verbose,
                    )
                else:

                    results[b0].append(this_date)

    return results


def clean(df):
    """
    pseudocode

    for each star
        for each timestamp:

            for each band:

                ask:
                    * is the data point more than 3sigma away from both its neighbors on each side?
                    * is its grade below 95% (or some other threshold)?
                    * do the other two bands haev data "within" 3 sigma of their neighbors?


    so there's probably a way in which I could code up some metric all in one pass for 
    the dataset which essentially assesses the first one -
    in other words, "compute how far in sigmas this data point is from its nearest 
    neighbors on either side."

    """
    pass
