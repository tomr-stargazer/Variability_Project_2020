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
