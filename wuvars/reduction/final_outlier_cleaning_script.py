"""
This script ~does~ the final outlier clean. 

It also helps me prototype the process.

"""

import matplotlib.pyplot as plt
import numpy as np
from final_outlier_clean import (calculate_diff,
                                 simple_lc_scatter_brokenaxes_grade)
from wuvars.data import photometry, quality_classes, spreadsheet


def mark_one_lightcurve():
    """
    This script creates one lightcurve that has its `problematic` data points marked
    with red circles.

    """

    # draw the lightcurve
    # (I am partial to the lightcurve with grade colors, because they are so salient
    #  to what I'm doing here.)

    dat = photometry.group_wserv_v2(photometry.load_wserv_v2(7))
    sid = 44508746107212
    date_offset = 56141

    sdat = dat.groups[dat.groups.keys["SOURCEID"] == sid]

    fig = simple_lc_scatter_brokenaxes_grade(
        dat, sid, date_offset=56141, xlims=[(-5, 251.0)]
    )

    # mark its problematic points with red circles

    # grab the J band ones
    # plot em on J

    j_outlier = np.max(sdat["JAPERMAG3"])
    date_outlier = sdat["MEANMJDOBS"][sdat["JAPERMAG3"] == j_outlier] - date_offset

    print(date_outlier, j_outlier)

    fig.ax_j.plot(
        [date_outlier],
        [j_outlier],
        ms=15,
        marker="o",
        markeredgecolor="r",
        markerfacecolor="None",
        zorder=100,
    )

    fig.canvas.draw()

    # (so, implicitly, we need a list of them)

    return fig


if __name__ == "__main__":

    fig = mark_one_lightcurve()
