"""
This script ~does~ the final outlier clean. 

It also helps me prototype the process.

"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.reduction.final_outlier_clean import (
    calculate_diffs, identify_outliers, simple_lc_scatter_brokenaxes_grade)


def mark_one_lightcurve(dat, sid, date_offset=56141, xlims=[(-5, 251.0)]):
    """
    This script creates one lightcurve that has its `problematic` data points marked
    with red circles.

    """

    # draw the lightcurve
    # (I am partial to the lightcurve with grade colors, because they are so salient
    #  to what I'm doing here.)

    sdat = dat.groups[dat.groups.keys["SOURCEID"] == sid]

    fig = simple_lc_scatter_brokenaxes_grade(
        dat, sid, date_offset=date_offset, xlims=xlims
    )

    # calculate this guy's problematic things

    diffs = calculate_diffs(dat, sid)

    # mark its problematic points with red circles

    # grab the J band ones
    # plot em on J

    outliers = identify_outliers(dat, sid, diffs, date_offset=date_offset)

    if False:
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
    else:

        j_outliers, h_outliers, k_outliers = outliers["J"], outliers["H"], outliers["K"]

        for j_date in j_outliers:

            j_date_offset = j_date - date_offset
            j_outlier = sdat["JAPERMAG3"][sdat["MEANMJDOBS"] == j_date]

            fig.ax_j.plot(
                [j_date_offset],
                [j_outlier],
                ms=15,
                marker="o",
                markeredgecolor="r",
                markerfacecolor="None",
                zorder=100,
            )

        for h_date in h_outliers:

            h_date_offset = h_date - date_offset
            h_outlier = sdat["HAPERMAG3"][sdat["MEANMJDOBS"] == h_date]

            fig.ax_h.plot(
                [h_date_offset],
                [h_outlier],
                ms=15,
                marker="o",
                markeredgecolor="r",
                markerfacecolor="None",
                zorder=100,
            )

        for k_date in k_outliers:

            k_date_offset = k_date - date_offset
            k_outlier = sdat["KAPERMAG3"][sdat["MEANMJDOBS"] == k_date]

            fig.ax_k.plot(
                [k_date_offset],
                [k_outlier],
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

    ngc_dat = photometry.group_wserv_v2(photometry.load_wserv_v2(7))
    new_data = photometry.group_wserv_v2(
        photometry.load_wserv_v2(7, suffix="_outlier_cleaned_152"), max_flags=257
    )

    ngc_date_offset = 56141
    ngc_xlims = [(-5, 251.0)]

    sid = 44508746107212
    sid2 = 44508746127611
    sid3 = 44508746107259
    sid4 = 44508746116310
    sid5 = 44508746098496

    sid_list = [sid, sid2, sid3, sid4, sid5]

    from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc

    # match_ngc.approved
    if True:
        for sidd in sid_list:

            fig = mark_one_lightcurve(
                ngc_dat, sidd, date_offset=ngc_date_offset, xlims=ngc_xlims
            )

            fig2 = mark_one_lightcurve(
                new_data, sidd, date_offset=ngc_date_offset, xlims=ngc_xlims
            )

    ic_match = match_ic()
    ic_dat = photometry.group_wserv_v2(photometry.load_wserv_v2(8))
    from wuvars.plotting.lightcurve import ic_date_offset, ic348_xlims

    ic_sids = ic_match.approved["SOURCEID"]

    if False:
        for i, sidd in enumerate(ic_sids):

            fig = mark_one_lightcurve(
                ic_dat, sidd, date_offset=ic_date_offset, xlims=ic348_xlims
            )

            if i > 10:
                break
