"""
I plan to put some light curve plotting code here.

"""

import numpy as np
import matplotlib.pyplot as plt


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

    ax_j.scatter(date, j, c=date, s=18, edgecolors='k', linewidths=0.5, **kwargs)
    ax_h.scatter(date, h, c=date, s=18, edgecolors='k', linewidths=0.5, **kwargs)
    ax_k.scatter(date, k, c=date, s=18, edgecolors='k', linewidths=0.5, **kwargs)

    ax_j.errorbar(date, j, yerr=j_e, fmt="None", ecolor="k", ms=2, elinewidth=0.5, zorder=-1)
    ax_h.errorbar(date, h, yerr=h_e, fmt="None", ecolor="k", ms=2, elinewidth=0.5, zorder=-1)
    ax_k.errorbar(date, k, yerr=k_e, fmt="None", ecolor="k", ms=2, elinewidth=0.5, zorder=-1)

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
        phase - 1, x, yerr=xerr, fmt=sym, mfc="0.7", mec="0.7", ecolor="0.7", ms=ms
    )
    ax.errorbar(
        phase + 1, x, yerr=xerr, fmt=sym, mfc="0.7", mec="0.7", ecolor="0.7", ms=ms
    )

    # Now plot our actual scattered dude
    if not hide:
        # errorbars in the background
        ax.errorbar(phase, x, yerr=xerr, fmt="None", ecolor="k", elinewidth=0.5, zorder=0)
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

    scatter_phase_core(ax_j, date, j, j_e, period, offset=offset, c=date, s=18, edgecolors='k', linewidths=0.5, ms=4, **kwargs)
    scatter_phase_core(ax_h, date, h, h_e, period, offset=offset, c=date, s=18, edgecolors='k', linewidths=0.5, ms=4, **kwargs)
    sc = scatter_phase_core(ax_k, date, k, k_e, period, offset=offset, c=date, s=18, edgecolors='k', linewidths=0.5, ms=4, **kwargs)

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

    cbar = plt.gcf().colorbar(sc, ax=ax_khk) # This should really be changed to the method
    cbar.set_label(f"MJD - {int(begin)}")

    cbar = plt.gcf().colorbar(sc, ax=ax_jhk) # This should really be changed to the method
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
