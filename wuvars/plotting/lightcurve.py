"""
I plan to put some light curve plotting code here.

"""

import matplotlib.pyplot as plt
import numpy as np
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec
from wuvars.plotting.lightcurve_helpers import orion_cmap, produce_xlims

onc_date_offset = 54034.0
ngc_date_offset = 56141
ic_date_offset = 56849
monr2_date_offset = 57374


def simple_lc(dg, sid):
    # dg: astropy table that has been grouped by SOURCEID

    dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

    # set up data

    date = dat["MEANMJDOBS"]

    j = dat["JAPERMAG3"]
    h = dat["HAPERMAG3"]
    k = dat["KAPERMAG3"]

    j_e = dat["JAPERMAG3ERR"]
    h_e = dat["HAPERMAG3ERR"]
    k_e = dat["KAPERMAG3ERR"]

    # set up plot

    fig = plt.figure(figsize=(10, 6), dpi=80, facecolor="w", edgecolor="k")

    bottom = 0.1
    height = 0.25
    left = 0.075
    width = 0.5

    ax_k = fig.add_axes((left, bottom, width, height))
    ax_h = fig.add_axes((left, bottom + 0.3, width, height), sharex=ax_k)
    ax_j = fig.add_axes((left, bottom + 0.6, width, height), sharex=ax_k)

    ax_jhk = fig.add_axes((0.65, bottom, 0.3, 0.375))
    ax_khk = fig.add_axes((0.65, bottom + 0.475, 0.3, 0.375))

    fig.ax_k = ax_k
    fig.ax_j = ax_j
    fig.ax_h = ax_h
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    ax_j.errorbar(date, j, yerr=j_e, fmt="bo", ecolor="k", ms=2, elinewidth=0.5)
    ax_h.errorbar(date, h, yerr=h_e, fmt="go", ecolor="k", ms=2, elinewidth=0.5)
    ax_k.errorbar(date, k, yerr=k_e, fmt="ro", ecolor="k", ms=2, elinewidth=0.5)

    ax_jhk.errorbar(
        h - k,
        j - h,
        xerr=(h_e ** 2 + k_e ** 2) ** 0.5,
        yerr=(h_e ** 2 + j_e ** 2) ** 0.5,
        fmt="k.",
        alpha=0.1,
    )

    ax_khk.errorbar(
        h - k, k, xerr=(h_e ** 2 + k_e ** 2) ** 0.5, yerr=k_e, fmt="k.", alpha=0.1
    )

    ax_j.invert_yaxis()
    ax_h.invert_yaxis()
    ax_k.invert_yaxis()
    ax_khk.invert_yaxis()

    ax_j.set_ylabel("J", {"rotation": "horizontal", "fontsize": "large"})
    ax_h.set_ylabel("H", {"rotation": "horizontal", "fontsize": "large"})
    ax_k.set_ylabel("K", {"rotation": "horizontal", "fontsize": "large"})

    ax_k.set_xlabel("MJD")

    ax_jhk.set_xlabel("H-K")
    ax_jhk.set_ylabel("J-H")  # , {'rotation':'horizontal'})
    ax_khk.set_xlabel("H-K")
    ax_khk.set_ylabel("K")  # , {'rotation':'horizontal'})

    return fig


def simple_lc_scatter(dg, sid, begin=0, **kwargs):
    # dg: astropy table that has been grouped by SOURCEID

    dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

    # set up data

    date = dat["MEANMJDOBS"] - begin

    j = dat["JAPERMAG3"]
    h = dat["HAPERMAG3"]
    k = dat["KAPERMAG3"]

    j_e = dat["JAPERMAG3ERR"]
    h_e = dat["HAPERMAG3ERR"]
    k_e = dat["KAPERMAG3ERR"]

    # set up plot

    fig = plt.figure(figsize=(10, 6), dpi=80, facecolor="w", edgecolor="k")

    bottom = 0.1
    height = 0.25
    left = 0.075
    width = 0.5

    ax_k = fig.add_axes((left, bottom, width, height))
    ax_h = fig.add_axes((left, bottom + 0.3, width, height), sharex=ax_k)
    ax_j = fig.add_axes((left, bottom + 0.6, width, height), sharex=ax_k)

    ax_jhk = fig.add_axes((0.65, bottom, 0.3, 0.375))
    ax_khk = fig.add_axes((0.65, bottom + 0.475, 0.3, 0.375))

    fig.ax_k = ax_k
    fig.ax_j = ax_j
    fig.ax_h = ax_h
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    ax_j.scatter(
        date,
        j,
        c=date,
        vmin=date.min(),
        vmax=date.max(),
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )
    ax_h.scatter(
        date,
        h,
        c=date,
        vmin=date.min(),
        vmax=date.max(),
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )
    ax_k.scatter(
        date,
        k,
        c=date,
        vmin=date.min(),
        vmax=date.max(),
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )

    ax_j.errorbar(
        date, j, yerr=j_e, fmt="None", ecolor="k", ms=2, elinewidth=0.5, zorder=-1
    )
    ax_h.errorbar(
        date, h, yerr=h_e, fmt="None", ecolor="k", ms=2, elinewidth=0.5, zorder=-1
    )
    ax_k.errorbar(
        date, k, yerr=k_e, fmt="None", ecolor="k", ms=2, elinewidth=0.5, zorder=-1
    )

    ax_jhk.errorbar(
        h - k,
        j - h,
        xerr=(h_e ** 2 + k_e ** 2) ** 0.5,
        yerr=(h_e ** 2 + j_e ** 2) ** 0.5,
        fmt="k.",
        alpha=0.1,
    )

    ax_khk.errorbar(
        h - k, k, xerr=(h_e ** 2 + k_e ** 2) ** 0.5, yerr=k_e, fmt="k.", alpha=0.1
    )

    ax_j.invert_yaxis()
    ax_h.invert_yaxis()
    ax_k.invert_yaxis()
    ax_khk.invert_yaxis()

    ax_j.set_ylabel("J", {"rotation": "horizontal", "fontsize": "large"})
    ax_h.set_ylabel("H", {"rotation": "horizontal", "fontsize": "large"})
    ax_k.set_ylabel("K", {"rotation": "horizontal", "fontsize": "large"})

    ax_k.set_xlabel(f"MJD - {int(begin)}")

    ax_jhk.set_xlabel("H-K")
    ax_jhk.set_ylabel("J-H")  # , {'rotation':'horizontal'})
    ax_khk.set_xlabel("H-K")
    ax_khk.set_ylabel("K")  # , {'rotation':'horizontal'})

    return fig


def simple_lc_brokenaxes(dg, sid, date_offset=None, pad=5, xlims=None, breaks=None):
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

    # set up plot
    # OUR brokenaxes changes will all go here
    # need:
    #   - how much MJD to subtract
    #   - where to put the breaks (+how many breaks)
    #   - ... that's p much it, right?

    fig = plt.figure(figsize=(10, 6), dpi=80, facecolor="w", edgecolor="k")

    gs0 = fig.add_gridspec(1, 2, width_ratios=[6, 3], wspace=0.2)
    gs_left = gs0[0].subgridspec(3, 1, hspace=0)
    gs_right = gs0[1].subgridspec(1, 2, width_ratios=(2, 0.1), wspace=0.15)
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
        print(xlims)

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

    ax_j.errorbar(date, j, yerr=j_e, fmt="bo", ecolor="k", ms=2, elinewidth=0.5)
    ax_h.errorbar(date, h, yerr=h_e, fmt="go", ecolor="k", ms=2, elinewidth=0.5)
    ax_k.errorbar(date, k, yerr=k_e, fmt="ro", ecolor="k", ms=2, elinewidth=0.5)

    ax_jhk.errorbar(
        h - k,
        j - h,
        xerr=(h_e ** 2 + k_e ** 2) ** 0.5,
        yerr=(h_e ** 2 + j_e ** 2) ** 0.5,
        fmt="k.",
        alpha=0.1,
    )

    ax_khk.errorbar(
        h - k, k, xerr=(h_e ** 2 + k_e ** 2) ** 0.5, yerr=k_e, fmt="k.", alpha=0.1
    )

    ax_j.invert_yaxis()
    ax_h.invert_yaxis()
    ax_k.invert_yaxis()
    ax_khk.invert_yaxis()

    ax_j.set_ylabel("J", labelpad=40, fontdict={"rotation": "horizontal"})
    ax_h.set_ylabel("H", labelpad=40, fontdict={"rotation": "horizontal"})
    ax_k.set_ylabel("K", labelpad=40, fontdict={"rotation": "horizontal"})

    # \u2212 is a proper minus sign (better than the hyphen character `-`)
    ax_k.set_xlabel(f"MJD \u2212 {date_offset}", labelpad=20)

    ax_jhk.set_xlabel("H-K")
    ax_jhk.set_ylabel("J-H")  # , {'rotation':'horizontal'})
    ax_khk.set_xlabel("H-K")
    ax_khk.set_ylabel("K")  # , {'rotation':'horizontal'})

    return fig


ic348_xlims = [(-5.0, 24.0), (74.0, 230.0), (370.0, 390.0)]


def ic348_simple_lc_brokenaxes(dg, sid, date_offset=ic_date_offset, xlims=ic348_xlims):

    return simple_lc_brokenaxes(dg, sid, date_offset=date_offset, xlims=xlims)


ngc1333_xlims = [(-5, 251.0)]


def ngc1333_simple_lc_brokenaxes(
    dg, sid, date_offset=ngc_date_offset, xlims=ngc1333_xlims
):

    return simple_lc_brokenaxes(dg, sid, date_offset=date_offset, xlims=xlims)


monr2_xlims = [(-5.0, 134.0), (291.0, 414.0)]


def monr2_simple_lc_brokenaxes(
    dg, sid, date_offset=monr2_date_offset, xlims=monr2_xlims
):

    return simple_lc_brokenaxes(dg, sid, date_offset=date_offset, xlims=xlims)


onc_xlims = [(-7.0, 185.0), (384.0, 410.0), (737.0, 756.0), (822.0, 902.0)]


def onc_simple_lc_brokenaxes(dg, sid, date_offset=onc_date_offset, xlims=onc_xlims):

    return simple_lc_brokenaxes(dg, sid, date_offset=date_offset, xlims=xlims)


def simple_lc_scatter_brokenaxes(
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
        c=date,
        vmin=date.min(),
        vmax=date.max(),
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )
    sc_h = ax_h.scatter(
        date,
        h,
        c=date,
        vmin=date.min(),
        vmax=date.max(),
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )
    sc_k = ax_k.scatter(
        date,
        k,
        c=date,
        vmin=date.min(),
        vmax=date.max(),
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

    ax_jhk.scatter(
        h - k,
        j - h,
        c=date,
        vmin=date.min(),
        vmax=date.max(),
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )

    ax_khk.scatter(
        h - k,
        k,
        c=date,
        vmin=date.min(),
        vmax=date.max(),
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

    return fig


def ic348_simple_lc_scatter_brokenaxes(
    dg, sid, date_offset=ic_date_offset, xlims=ic348_xlims, **kwargs
):

    return simple_lc_scatter_brokenaxes(
        dg, sid, date_offset=date_offset, xlims=xlims, **kwargs
    )


def ngc1333_simple_lc_scatter_brokenaxes(
    dg, sid, date_offset=ngc_date_offset, xlims=ngc1333_xlims, **kwargs
):

    return simple_lc_scatter_brokenaxes(
        dg, sid, date_offset=date_offset, xlims=xlims, **kwargs
    )


def monr2_simple_lc_scatter_brokenaxes(
    dg, sid, date_offset=monr2_date_offset, xlims=monr2_xlims, **kwargs
):

    return simple_lc_scatter_brokenaxes(
        dg, sid, date_offset=date_offset, xlims=xlims, **kwargs
    )


# Note: This one is special, because it has a default cmap assigned
def onc_simple_lc_scatter_brokenaxes(
    dg, sid, date_offset=onc_date_offset, xlims=onc_xlims, cmap=orion_cmap, **kwargs
):

    return simple_lc_scatter_brokenaxes(
        dg, sid, date_offset=date_offset, xlims=xlims, cmap=cmap, **kwargs
    )


# 'hide' parameter is hackish.
def plot_phase_core(
    ax, t, x, xerr, period, offset=0, sym="o", color="k", ms=4, hide=False
):
    """
    Plots a pretty period-folded lightcurve on a given axes object.

    Doesn't assume anything about your data (e.g., that it's in magnitudes)

    Parameters
    ----------
    ax : plt.Axes
    t, x, xerr : array_like
    period : float
    offset : float, optional
        How much to shift the phase by. Default is zero.
    sym : str, optional
        Default 'o'. (circles)
    color : str, optional
        Default 'k'. (black)
    ms : float
        Default 6.

    Returns
    -------
    period : float
        The input period.

    """

    phase = ((t % period) / period + offset) % 1.0

    if not hide:
        ax.errorbar(phase, x, yerr=xerr, fmt=color + sym, ms=ms)
    ax.errorbar(
        phase - 1, x, yerr=xerr, fmt=sym, mfc="0.7", mec="0.7", ecolor="0.7", ms=ms
    )
    ax.errorbar(
        phase + 1, x, yerr=xerr, fmt=sym, mfc="0.7", mec="0.7", ecolor="0.7", ms=ms
    )

    ax.set_xticks([0, 0.5, 1])
    ax.set_xticks(np.arange(-0.5, 1.5, 0.1), minor=True)

    ax.set_xlim(-0.25, 1.25)

    return period


def scatter_phase_core(
    ax, t, x, xerr, period, offset=0, ms=6, sym="o", hide=False, **kwargs
):
    """ 
    Scatter-plots a period-folded lightcurve on a given axes object.

    Also plots errorbars underneath.
    Doesn't assume anything about your data (e.g., that it's in magnitudes)
    
    Parameters
    ----------
    ax : plt.Axes
        An Axes object to plot onto.
    t, x, xerr : array_like
        The time, value, and uncertainty arrays to plot.
    period : float
        The period to fold the curve by.
    offset : float, optional
        How much to shift the phase by. Default is zero.
    sym : str, optional
        Symbol for gray dudes on the sides. Default 'o'. (circles)
    color : str, optional
        Color of the errorbars. Default 'k'. (black)
    ms : float
        Default 6.
    **kwargs : keyword arguments for plt.scatter()
        Suggested kwargs: `c` for color-array, `cmap` for colormap choice,
        `vmin` and `vmax` for the limits of colormap. Also possible:
        s (size), other things too I guess.
        
    Returns
    -------
    period : float
        The input period, unchanged.
    
    """

    # Calculate the "phase" variable as a function of time and the period.
    phase = ((t % period) / period + offset) % 1.0

    # Plot the grayed-out guys on the left and right:
    ax.errorbar(
        phase - 1,
        x,
        yerr=xerr,
        fmt=sym,
        mfc="0.7",
        mec="0.7",
        ecolor="0.7",
        ms=ms,
        alpha=0.5,
    )
    ax.errorbar(
        phase + 1,
        x,
        yerr=xerr,
        fmt=sym,
        mfc="0.7",
        mec="0.7",
        ecolor="0.7",
        ms=ms,
        alpha=0.5,
    )

    # Now plot our actual scattered dude
    if not hide:
        # errorbars in the background
        ax.errorbar(
            phase,
            x,
            yerr=xerr,
            fmt="None",
            ecolor="k",
            elinewidth=0.5,
            zorder=0,
            alpha=0.5,
        )
        # scatter in the foreground
        sc = ax.scatter(phase, x, zorder=100, **kwargs)

    ax.set_xticks([0, 0.5, 1])
    ax.set_xticks(np.arange(-0.5, 1.5, 0.1), minor=True)

    ax.set_xlim(-0.1, 1.1)

    return sc


def simple_phased_lc(dg, sid, period, offset=0):
    # dg: astropy table that has been grouped by SOURCEID

    dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

    # set up data

    date = dat["MEANMJDOBS"]

    j = dat["JAPERMAG3"]
    h = dat["HAPERMAG3"]
    k = dat["KAPERMAG3"]

    j_e = dat["JAPERMAG3ERR"]
    h_e = dat["HAPERMAG3ERR"]
    k_e = dat["KAPERMAG3ERR"]

    # set up plot

    fig = plt.figure(figsize=(10, 6), dpi=80, facecolor="w", edgecolor="k")

    bottom = 0.1
    height = 0.25
    left = 0.075
    width = 0.5

    ax_k = fig.add_axes((left, bottom, width, height))
    ax_h = fig.add_axes((left, bottom + 0.3, width, height), sharex=ax_k)
    ax_j = fig.add_axes((left, bottom + 0.6, width, height), sharex=ax_k)

    ax_jhk = fig.add_axes((0.65, bottom, 0.3, 0.375))
    ax_khk = fig.add_axes((0.65, bottom + 0.475, 0.3, 0.375))

    fig.ax_k = ax_k
    fig.ax_j = ax_j
    fig.ax_h = ax_h
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    plot_phase_core(ax_j, date, j, j_e, period, offset=offset, color="b")
    plot_phase_core(ax_h, date, h, h_e, period, offset=offset, color="g")
    plot_phase_core(ax_k, date, k, k_e, period, offset=offset, color="r")

    ax_jhk.errorbar(
        h - k,
        j - h,
        xerr=(h_e ** 2 + k_e ** 2) ** 0.5,
        yerr=(h_e ** 2 + j_e ** 2) ** 0.5,
        fmt="k.",
        alpha=0.1,
    )

    ax_khk.errorbar(
        h - k, k, xerr=(h_e ** 2 + k_e ** 2) ** 0.5, yerr=k_e, fmt="k.", alpha=0.1
    )

    ax_j.invert_yaxis()
    ax_h.invert_yaxis()
    ax_k.invert_yaxis()
    ax_khk.invert_yaxis()

    ax_j.set_ylabel("J", {"rotation": "horizontal", "fontsize": "large"})
    ax_h.set_ylabel("H", {"rotation": "horizontal", "fontsize": "large"})
    ax_k.set_ylabel("K", {"rotation": "horizontal", "fontsize": "large"})

    ax_k.set_xlabel(f"Phase (Period={period:.2f}d)")

    ax_jhk.set_xlabel("H-K")
    ax_jhk.set_ylabel("J-H")  # , {'rotation':'horizontal'})
    ax_khk.set_xlabel("H-K")
    ax_khk.set_ylabel("K")  # , {'rotation':'horizontal'})

    return fig


def simple_phased_lc_scatter(dg, sid, period, offset=0, begin=0, **kwargs):
    # dg: astropy table that has been grouped by SOURCEID

    dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

    # set up data

    date = dat["MEANMJDOBS"] - begin

    j = dat["JAPERMAG3"]
    h = dat["HAPERMAG3"]
    k = dat["KAPERMAG3"]

    j_e = dat["JAPERMAG3ERR"]
    h_e = dat["HAPERMAG3ERR"]
    k_e = dat["KAPERMAG3ERR"]

    # set up plot

    fig = plt.figure(figsize=(10, 6), dpi=80, facecolor="w", edgecolor="k")

    bottom = 0.1
    height = 0.25
    left = 0.075
    width = 0.5

    ax_k = fig.add_axes((left, bottom, width, height))
    ax_h = fig.add_axes((left, bottom + 0.3, width, height), sharex=ax_k)
    ax_j = fig.add_axes((left, bottom + 0.6, width, height), sharex=ax_k)

    ax_jhk = fig.add_axes((0.65, bottom, 0.3, 0.375))
    ax_khk = fig.add_axes((0.65, bottom + 0.475, 0.3, 0.375))

    fig.ax_k = ax_k
    fig.ax_j = ax_j
    fig.ax_h = ax_h
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    scatter_phase_core(
        ax_j,
        date,
        j,
        j_e,
        period,
        offset=offset,
        c=date,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        ms=4,
        **kwargs,
    )
    scatter_phase_core(
        ax_h,
        date,
        h,
        h_e,
        period,
        offset=offset,
        c=date,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        ms=4,
        **kwargs,
    )
    sc = scatter_phase_core(
        ax_k,
        date,
        k,
        k_e,
        period,
        offset=offset,
        c=date,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        ms=4,
        **kwargs,
    )

    ax_jhk.errorbar(
        h - k,
        j - h,
        xerr=(h_e ** 2 + k_e ** 2) ** 0.5,
        yerr=(h_e ** 2 + j_e ** 2) ** 0.5,
        fmt="k.",
        alpha=0.1,
    )

    ax_khk.errorbar(
        h - k, k, xerr=(h_e ** 2 + k_e ** 2) ** 0.5, yerr=k_e, fmt="k.", alpha=0.1
    )

    ax_j.invert_yaxis()
    ax_h.invert_yaxis()
    ax_k.invert_yaxis()
    ax_khk.invert_yaxis()

    ax_j.set_ylabel("J", {"rotation": "horizontal", "fontsize": "large"})
    ax_h.set_ylabel("H", {"rotation": "horizontal", "fontsize": "large"})
    ax_k.set_ylabel("K", {"rotation": "horizontal", "fontsize": "large"})

    ax_k.set_xlabel(f"Phase (Period={period:.2f}d)")

    ax_jhk.set_xlabel("H-K")
    ax_jhk.set_ylabel("J-H")  # , {'rotation':'horizontal'})
    ax_khk.set_xlabel("H-K")
    ax_khk.set_ylabel("K")  # , {'rotation':'horizontal'})

    cbar = plt.gcf().colorbar(
        sc, ax=ax_khk
    )  # This should really be changed to the method
    cbar.set_label(f"MJD - {int(begin)}")

    cbar = plt.gcf().colorbar(
        sc, ax=ax_jhk
    )  # This should really be changed to the method
    cbar.set_label(f"MJD - {int(begin)}")

    return fig


# from plot4 - we'll transform this
def phase_axes_with_info(
    stardata, band, period, axes, colorscale, cmap, vmin, vmax, offset=0, **kwargs
):

    columns = stardata.get_columns(band, max_flag=0)
    columns_info = stardata.get_columns(band, min_flag=1, max_flag=256)

    date = np.copy(columns["date"])
    date_info = np.copy(columns_info["date"])

    phase = ((date % period) / period + offset) % 1
    phase_info = ((date_info % period) / period + offset) % 1

    if len(columns["date"]) > 0:
        # plot the greyed-out versions on left and right
        axes.errorbar(
            phase - 1,
            columns["mag"],
            yerr=columns["err"],
            fmt="o",
            mfc="0.7",
            mec="0.7",
            ecolor="0.7",
            ms=6,
            zorder=-5,
        )
        axes.errorbar(
            phase + 1,
            columns["mag"],
            yerr=columns["err"],
            fmt="o",
            mfc="0.7",
            mec="0.7",
            ecolor="0.7",
            ms=6,
            zorder=-5,
        )

        # First, plot the errorbars, with no markers, in the background:
        axes.errorbar(
            phase,
            columns["mag"],
            marker=None,
            yerr=columns["err"],
            fmt="none",
            ecolor="k",
            zorder=0,
        )
        # Next, scatter the points themselves, colored re:colorscale :
        axes.scatter(
            phase,
            columns["mag"],
            cmap=cmap,
            c=columns[colorscale],
            vmin=vmin,
            vmax=vmax,
            zorder=100,
            **kwargs,
        )

    if len(columns_info["date"]) > 0:
        # plot the greyed-out versions on left and right
        axes.errorbar(
            phase_info - 1,
            columns_info["mag"],
            yerr=columns_info["err"],
            fmt="d",
            mfc="0.7",
            mec="0.7",
            ecolor="0.7",
            ms=6,
            zorder=-5,
        )
        axes.errorbar(
            phase_info + 1,
            columns_info["mag"],
            yerr=columns_info["err"],
            fmt="d",
            mfc="0.7",
            mec="0.7",
            ecolor="0.7",
            ms=6,
            zorder=-5,
        )

        # First, plot the errorbars, with no markers, in the background:
        axes.errorbar(
            phase_info,
            columns_info["mag"],
            yerr=columns_info["err"],
            marker=None,
            fmt="none",
            ecolor="k",
            zorder=0,
        )
        # Next, scatter the points themselves, colored re:colorscale :
        axes.scatter(
            phase_info,
            columns_info["mag"],
            marker="d",
            c=columns_info[colorscale],
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            zorder=100,
            **kwargs,
        )

    # Finally, flip it (magnitudes are backwards).
    axes.invert_yaxis()

    axes.set_xticks([0, 0.5, 1])
    axes.set_xticks(np.arange(-0.5, 1.5, 0.1), minor=True)

    axes.set_xlim(-0.25, 1.25)

    axes.get_figure().canvas.draw()


def simple_phased_lc_scatter_gridspec(
    dg, sid, period, date_offset=None, phase_offset=0, color_by="date", **kwargs
):
    # dg: astropy table that has been grouped by SOURCEID

    dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

    # set up data

    if date_offset is None:
        date_offset = np.floor(np.min(dat["MEANMJDOBS"]))
    date = dat["MEANMJDOBS"] - date_offset
    phase = ((date % period) / period + phase_offset) % 1.0

    if color_by == "date":
        color_array = date
    elif color_by == "phase":
        color_array = phase
    else:
        raise ValueError("`color_by` option must be 'date' or 'phase'")

    j = dat["JAPERMAG3"]
    h = dat["HAPERMAG3"]
    k = dat["KAPERMAG3"]

    j_e = dat["JAPERMAG3ERR"]
    h_e = dat["HAPERMAG3ERR"]
    k_e = dat["KAPERMAG3ERR"]

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

    ax_j = fig.add_subplot(gs_left[0, 0])
    ax_h = fig.add_subplot(gs_left[1, 0])
    ax_k = fig.add_subplot(gs_left[2, 0])

    ax_khk = fig.add_subplot(gs01[0, 0])
    ax_jhk = fig.add_subplot(gs01[1, 0])

    # ax_cbar = fig.add_subplot(gs_right[1])

    ax_j.tick_params(labelbottom=False)
    ax_h.tick_params(labelbottom=False)

    fig.ax_k = ax_k
    fig.ax_j = ax_j
    fig.ax_h = ax_h
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    scatter_phase_core(
        ax_j,
        date,
        j,
        j_e,
        period,
        offset=phase_offset,
        c=color_array,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        ms=4,
        **kwargs,
    )
    scatter_phase_core(
        ax_h,
        date,
        h,
        h_e,
        period,
        offset=phase_offset,
        c=color_array,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        ms=4,
        **kwargs,
    )
    sc = scatter_phase_core(
        ax_k,
        date,
        k,
        k_e,
        period,
        offset=phase_offset,
        c=color_array,
        s=18,
        edgecolors="k",
        linewidths=0.5,
        ms=4,
        **kwargs,
    )

    ax_jhk.scatter(
        h - k,
        j - h,
        c=color_array,
        vmin=color_array.min(),
        vmax=color_array.max(),
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )

    ax_khk.scatter(
        h - k,
        k,
        c=color_array,
        vmin=color_array.min(),
        vmax=color_array.max(),
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

    ax_j.invert_yaxis()
    ax_h.invert_yaxis()
    ax_k.invert_yaxis()
    ax_khk.invert_yaxis()

    ax_j.set_ylabel("J", {"rotation": "horizontal", "fontsize": "large"})
    ax_h.set_ylabel("H", {"rotation": "horizontal", "fontsize": "large"})
    ax_k.set_ylabel("K", {"rotation": "horizontal", "fontsize": "large"})

    ax_k.set_xlabel(f"Phase (Period={period:.2f}d)")

    ax_jhk.set_xlabel("H-K")
    ax_jhk.set_ylabel("J-H")  # , {'rotation':'horizontal'})
    ax_khk.set_xlabel("H-K")
    ax_khk.set_ylabel("K")  # , {'rotation':'horizontal'})

    return fig


# Trying something!
def cody_singleax_lc_brokenaxes(
    dg, sid, date_offset=None, pad=5, xlims=None, breaks=None
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

    # set up plot
    # OUR brokenaxes changes will all go here
    # need:
    #   - how much MJD to subtract
    #   - where to put the breaks (+how many breaks)
    #   - ... that's p much it, right?

    fig = plt.figure(figsize=(5, 1), dpi=240, facecolor="w", edgecolor="k")
    bax_kwargs = dict(despine=False, d=0.0075, tilt=60, wspace=0.04)

    if xlims is None:
        if breaks is None:
            # baked-in default since I was prototyping on IC348; possibly a bad default
            breaks = [50, 350]
        xlims = produce_xlims(date, breaks=breaks, pad=pad)
        print(xlims)

    ax_j = brokenaxes(xlims=xlims, fig=fig, **bax_kwargs)

    # Make the broken zones gray. (Uses some details from  )
    ax_j.big_ax.set_zorder(-100)
    ax_j.big_ax.set_facecolor("0.9")

    fig.ax_j = ax_j

    offset = k.mean() - j.mean()
    # ax_j.errorbar(date, h, yerr=h_e, fmt="go", ecolor="k", ms=2, elinewidth=0.5)
    ax_j.errorbar(date, j, yerr=j_e, fmt="k.", ecolor="k", ms=2, elinewidth=0.25)
    ax_j.errorbar(
        date,
        k - offset,
        yerr=k_e,
        fmt=".",
        color="0.75",
        ecolor="0.75",
        ms=5,
        elinewidth=0.5,
        zorder=-1,
    )
    ax_j.invert_yaxis()

    ax_j.set_ylabel("J", labelpad=40, fontdict={"rotation": "horizontal"})

    # \u2212 is a proper minus sign (better than the hyphen character `-`)
    ax_j.set_xlabel(f"MJD \u2212 {date_offset}", labelpad=20)

    return fig


def trim_bottom_labels(ax, remove_yticklabel=False, threshold=0.1):
    ax.tick_params(labelbottom=False)

    if remove_yticklabel:
        # Get bottom y tick label
        labels = ax.get_yticklabels()
        try:
            labels[0].set_visible(False)
        except AttributeError:
            labels[0][0].set_visible(False)
    # if labels:
    #     bottom_label = labels[0][0] # not sure why it's twice

    #     # Get position in axis coordinates
    #     pos_data = bottom_label.get_position()  # (x, y) in data units

    #     try
    #     pos_pix = ax.transData.transform(pos_data)
    #     pos_axes = ax.transAxes.inverted().transform(pos_pix)

    #     # Check vertical coordinate (axis coords)
    #     if pos_axes[1] <= threshold:
    #         bottom_label.set_visible(False)

    # ax.get_yticklabels()[0].set_visible(False)

    return None


def eightpanel_lc(
    dg,
    sid,
    period,
    date_offset=None,
    phase_offset=0,
    color_by="date",
    pad=5,
    xlims=None,
    breaks=None,
    detrended_dg=None,
    **kwargs,
):
    """
    This is a joint 'straight-lc' and 'phase-folded-lc' in which the latter has detrending applied.

    This is intended for making a publication FigureSet.

    Borrows code largely from the following places:
    - simple_phased_lc_scatter_gridspec
    - simple_lc_scatter_brokenaxes

    """

    dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

    # data setup part 1

    if date_offset is None:
        date_offset = np.floor(np.min(dat["MEANMJDOBS"]))
    date = dat["MEANMJDOBS"] - date_offset
    phase = ((date % period) / period + phase_offset) % 1.0

    if color_by == "date":
        color_array = date
    elif color_by == "phase":
        color_array = phase
    else:
        raise ValueError("`color_by` option must be 'date' or 'phase'")

    if detrended_dg is None:
        detrended_dg = dat

    # fig setup

    fig = plt.figure(figsize=(20, 6), dpi=80, facecolor="w", edgecolor="k")

    gs0 = fig.add_gridspec(1, 3, width_ratios=[6, 6, 3], wspace=0.2)

    gs_left = gs0[0].subgridspec(3, 1, hspace=0)
    gs_center = gs0[1].subgridspec(3, 1, hspace=0)
    gs_right = gs0[2].subgridspec(1, 2, width_ratios=(2, 0.075), wspace=0.15)

    gs01 = gs_right[0].subgridspec(2, 1, hspace=0)

    # ax_j = fig.add_subplot(gs_left[0, 0])
    # ax_h = fig.add_subplot(gs_left[1, 0])
    # ax_k = fig.add_subplot(gs_left[2, 0])

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

    ax_j = brokenaxes(xlims=xlims, subplot_spec=gs_left[0, 0], **bax_kwargs)
    ax_h = brokenaxes(xlims=xlims, subplot_spec=gs_left[1, 0], **bax_kwargs)
    ax_k = brokenaxes(xlims=xlims, subplot_spec=gs_left[2, 0], **bax_kwargs)

    # Make the broken zones gray. (Uses some details from ...? )
    ax_j.big_ax.set_zorder(-100)
    ax_h.big_ax.set_zorder(-100)
    ax_k.big_ax.set_zorder(-100)
    ax_j.big_ax.set_facecolor("0.9")
    ax_h.big_ax.set_facecolor("0.9")
    ax_k.big_ax.set_facecolor("0.9")

    ax_j_phase = fig.add_subplot(gs_center[0, 0])
    ax_h_phase = fig.add_subplot(gs_center[1, 0])
    ax_k_phase = fig.add_subplot(gs_center[2, 0])

    ax_khk = fig.add_subplot(gs01[0, 0])
    ax_jhk = fig.add_subplot(gs01[1, 0])

    # data setup part 2

    bands = ["J", "H", "K"]

    mags = {}
    errs = {}
    mags_detrended = {}

    lc_axes = {"J": ax_j, "H": ax_h, "K": ax_k}
    phase_axes = {"J": ax_j_phase, "H": ax_h_phase, "K": ax_k_phase}

    for band in bands:
        mags[band] = dat[f"{band}APERMAG3"]
        errs[band] = dat[f"{band}APERMAG3ERR"]
        mags_detrended[band] = detrended_dg[f"{band}APERMAG3"]

    # this section is to start doing the plotting work

    # do the linear lc part (uses brokenaxes)

    for band in bands:

        lc_axes[band].scatter(
            date,
            mags[band],
            c=color_array,
            # vmin=date.min(),
            # vmax=date.max(),
            s=18,
            edgecolors="k",
            linewidths=0.5,
            **kwargs,
        )

        lc_axes[band].errorbar(
            date,
            mags[band],
            yerr=errs[band],
            fmt="None",
            ecolor="k",
            ms=2,
            elinewidth=0.5,
            zorder=-1,
            alpha=0.5,
        )

    # do the folded part

    for band in bands:

        scatter_phase_core(
            phase_axes[band],
            date,
            mags_detrended[band],
            errs[band],
            period,
            offset=phase_offset,
            c=color_array,
            s=18,
            edgecolors="k",
            linewidths=0.5,
            ms=4,
            **kwargs,
        )

    # do the colors

    j = mags["J"]
    h = mags["H"]
    k = mags["K"]

    j_e = errs["J"]
    h_e = errs["H"]
    k_e = errs["K"]

    ax_jhk.scatter(
        h - k,
        j - h,
        c=color_array,
        vmin=date.min(),
        vmax=date.max(),
        s=18,
        edgecolors="k",
        linewidths=0.5,
        **kwargs,
    )

    ax_khk.scatter(
        h - k,
        k,
        c=color_array,
        vmin=date.min(),
        vmax=date.max(),
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

    # final adjustments

    ax_j.invert_yaxis()
    ax_h.invert_yaxis()
    ax_k.invert_yaxis()

    ax_j_phase.invert_yaxis()
    ax_h_phase.invert_yaxis()
    ax_k_phase.invert_yaxis()

    ax_khk.invert_yaxis()

    ax_j.set_ylabel("J", labelpad=40, fontdict={"rotation": "horizontal"})
    ax_h.set_ylabel("H", labelpad=40, fontdict={"rotation": "horizontal"})
    ax_k.set_ylabel("K", labelpad=40, fontdict={"rotation": "horizontal"})

    ax_k.set_xlabel(f"MJD \u2212 {date_offset}", labelpad=20)

    ax_k_phase.set_xlabel(f"Phase (Period={period:.2f}d)")

    ax_jhk.set_ylabel("J-H")  # , {'rotation':'horizontal'})
    ax_jhk.set_xlabel("H-K")
    ax_khk.set_ylabel("K")  # , {'rotation':'horizontal'})

    trim_bottom_labels(ax_j, remove_yticklabel=True)
    trim_bottom_labels(ax_h, remove_yticklabel=True)

    trim_bottom_labels(ax_j_phase, remove_yticklabel=True)
    trim_bottom_labels(ax_h_phase, remove_yticklabel=True)

    trim_bottom_labels(ax_khk, remove_yticklabel=True)

    fig.ax_j = ax_j
    fig.ax_h = ax_h
    fig.ax_k = ax_k

    fig.ax_j_phase = ax_j_phase
    fig.ax_h_phase = ax_h_phase
    fig.ax_k_phase = ax_k_phase

    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    return fig
    

def ic348_eightpanel_lc(
    *args, date_offset=ic_date_offset, xlims=ic348_xlims, **kwargs
):

    return eightpanel_lc(
        *args, date_offset=date_offset, xlims=xlims, **kwargs
    )


def ngc1333_eightpanel_lc(
    *args, date_offset=ngc_date_offset, xlims=ngc1333_xlims, **kwargs
):

    return eightpanel_lc(
        *args, date_offset=date_offset, xlims=xlims, **kwargs
    )
