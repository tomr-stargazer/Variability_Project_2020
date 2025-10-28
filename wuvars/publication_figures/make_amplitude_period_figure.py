"""
Placeholder script for making some kind of amplitude-period figure.

Okay. Tom has decided to make a 3x3 panel figure which incorporates histograms.
This may be too much detail. But. We're gonna try.

Gonna steal a lot from run_v4_period_results

"""

import os
import pdb

import astropy.table
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.detrending import poly_detrend
from wuvars.analysis.load_periodics_v4 import (ic_periods, ngc_periods,
                                               select_periodic_variables_v4)
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.data.preferred_photometry import (photometry_wserv7,
                                              photometry_wserv8)
from wuvars.plotting.lightcurve import (eightpanel_lc, eightpanel_lc_v2,
                                        fivepanel_lc, ic348_eightpanel_lc,
                                        ic348_eightpanel_lc_v2,
                                        ic348_fivepanel_lc, ic_date_offset,
                                        ngc1333_eightpanel_lc,
                                        ngc1333_eightpanel_lc_v2,
                                        ngc1333_fivepanel_lc, ngc_date_offset)
from wuvars.publication_figures.make_figures import (figure_export_path,
                                                     figureset_export_path)

ngc_match = match_ngc()
ic_match = match_ic()


region_keys = ["ic", "ngc", "both"]
fullname_dict = {"ngc": "NGC 1333", "ic": "IC 348", 'both': "Combined"}
wserv_dict = {"ngc": 7, "ic": 8}

spread_dict = {"ngc": spreadsheet.load_wserv_v2(7), "ic": spreadsheet.load_wserv_v2(8)}
match_dict = {"ngc": ngc_match, "ic": ic_match}
per_dict = {"ngc": ngc_periods, "ic": ic_periods}

if __name__ == "__main__":

    # def make_the_figure():

    with mpl.rc_context(
        {
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "font.size": 6,
        }
    ):
        fig = plt.figure(figsize=(7.25, 7.5), dpi=150)

        gs = fig.add_gridspec(3, 3)
        # , wspace=0.2)

        # gs_left = gs0[0].subgridspec(3, 1, hspace=0)
        # gs_center = gs0[1].subgridspec(3, 1, hspace=0)
        # gs_right = gs0[2].subgridspec(1, 2, width_ratios=(2, 0.075), wspace=0.15)

        for i, region_key in enumerate(region_keys):

            # each region is a COLUMN so we're gonna

            if region_key == "both":
                periods = pd.concat([per_dict["ngc"], per_dict["ic"]], axis=0)

                # ngc_match_approved_copy = ngc_match.approved.copy()
                # ic_match_approved_copy = ic_match.approved.copy()

                # keep_columns = ['IRexc', 'SpT']

                # ngc_match_approved_copy.keep_columns(keep_columns)
                # ic_match_approved_copy.keep_columns(keep_columns)
                # approved = astropy.table.vstack([ngc_match_approved_copy, ic_match_approved_copy])
                approved = astropy.table.vstack([ngc_match.approved, ic_match.approved])
                spread = pd.concat([spread_dict["ngc"], spread_dict["ic"]], axis=0)

            else:
                periods = per_dict[region_key]
                approved = match_dict[region_key].approved
                spread = spread_dict[region_key]

            # let's divide up by infrared excess
            ir_exc = approved["IRexc"] == "yes"
            no_ir_exc = approved["IRexc"] == "no"
            neither_exc = (
                (approved["IRexc"] != "yes") & (approved["IRexc"] != "no")
            ).mask
            # print(
            #     "Neither:",
            #     len(approved)
            #     - np.sum(approved["IRexc"] == "yes")
            #     - np.sum(approved["IRexc"] == "no"),
            #     np.sum(neither_exc)
            #     # np.sum((approved["IRexc"] != "yes") & (approved["IRexc"] != "no"))
            # )

            per_nir_50 = np.nanmedian(periods["Period"][no_ir_exc])
            per_nir_84 = np.nanpercentile(periods["Period"][no_ir_exc], 84)
            per_nir_16 = np.nanpercentile(periods["Period"][no_ir_exc], 16)

            per_ir_50 = np.nanmedian(periods["Period"][ir_exc])
            per_ir_84 = np.nanpercentile(periods["Period"][ir_exc], 84)
            per_ir_16 = np.nanpercentile(periods["Period"][ir_exc], 16)

            amp_nir_50 = np.nanmedian(periods["Amp"][no_ir_exc])
            amp_nir_84 = np.nanpercentile(periods["Amp"][no_ir_exc], 84)
            amp_nir_16 = np.nanpercentile(periods["Amp"][no_ir_exc], 16)

            amp_ir_50 = np.nanmedian(periods["Amp"][ir_exc])
            amp_ir_84 = np.nanpercentile(periods["Amp"][ir_exc], 84)
            amp_ir_16 = np.nanpercentile(periods["Amp"][ir_exc], 16)

            #
            gs_cornery = gs[0, i].subgridspec(
                2, 2, width_ratios=(1, 0.4), height_ratios=(0.4, 1), hspace=0, wspace=0
            )

            ax_invisible = fig.add_subplot(gs[0, i])
            ax_invisible.set_axis_off()
            ax_invisible.set_title(fullname_dict[region_key])

            ax_amp_v_per = fig.add_subplot(gs_cornery[1, 0], zorder=100)
            ax_amp_hist = fig.add_subplot(gs_cornery[1, 1], sharey=ax_amp_v_per)
            ax_amp_hist.tick_params(labelleft=False)

            ax_per_hist = fig.add_subplot(gs_cornery[0, 0], sharex=ax_amp_v_per)
            ax_per_hist.tick_params(labelbottom=False)

            ax_amp_v_per.plot(
                periods["Period"][ir_exc],
                periods["Amp"][ir_exc],
                "rs",
                ms=1.5,
                label="IR excess",
            )
            ax_amp_v_per.plot(
                periods["Period"][no_ir_exc],
                periods["Amp"][no_ir_exc],
                "kD",
                ms=1.25,
                label="no IR excess",
            )
            # ax_amp_v_per.plot(
            #     periods["Period"][neither_exc],
            #     periods["Amp"][neither_exc],
            #     "o",
            #     color='0.5',
            #     ms=1.5,
            #     label="unknown IR excess",
            # )
            ax_amp_v_per.loglog()
            ax_amp_v_per.xaxis.set_major_formatter(mticker.ScalarFormatter())
            # ax_amp_v_per.yaxis.set_major_formatter(mticker.ScalarFormatter())

            ax_amp_v_per.set_xlabel("Period (d)")
            if i == 0:
                ax_amp_v_per.set_ylabel("Amplitude (mag)")

            min_per_bin = np.log10(np.nanmin(periods["Period"]))
            max_per_bin = np.log10(np.nanmax(periods["Period"]))
            # print("Log Per bins: ", min_per_bin, max_per_bin)

            min_amp_bin = np.log10(np.nanmin(periods["Amp"]))
            max_amp_bin = np.log10(np.nanmax(periods["Amp"]))
            # print("Log amp bins: ", min_amp_bin, max_amp_bin)

            min_per_bin = min(-0.11314264643119235, -0.2091950052836405) * 0.99
            max_per_bin = max(1.744725685332582, 1.2013475449906188) * 1.01

            min_amp_bin = min(-2.3159896476462887, -2.091253337507326) * 0.99
            max_amp_bin = max(-0.6408834417764375, -0.9407821623586071) * 1.01

            ax_per_hist.hist(
                periods["Period"][ir_exc],
                bins=np.logspace(min_per_bin, max_per_bin, 10),
                edgecolor="r",
                facecolor="None",
                hatch="///",
                histtype="stepfilled",
                label="IR excess",
            )
            ax_per_hist.hist(
                periods["Period"][~ir_exc],
                bins=np.logspace(min_per_bin, max_per_bin, 10),
                edgecolor="k",
                linestyle="--",
                facecolor="None",
                hatch="..",
                histtype="stepfilled",
                label="no IR excess",
            )

            ax_amp_hist.hist(
                periods["Amp"][ir_exc],
                bins=np.logspace(min_amp_bin, max_amp_bin, 10),
                orientation="horizontal",
                edgecolor="r",
                facecolor="None",
                hatch="///",
                histtype="stepfilled",
                label="IR excess",
            )

            ax_amp_hist.hist(
                periods["Amp"][~ir_exc],
                bins=np.logspace(min_amp_bin, max_amp_bin, 10),
                orientation="horizontal",
                edgecolor="k",
                linestyle="--",
                facecolor="None",
                hatch="..",
                histtype="stepfilled",
                label="no IR excess",
            )

            # ax_per_hist.axvline(per_ir_50, color="r", linestyle="--")
            # ax_per_hist.axvline(per_nir_50, color="k", linestyle=":")
            # ax_amp_hist.axhline(amp_ir_50, color="r", linestyle="--")
            # ax_amp_hist.axhline(amp_nir_50, color="k", linestyle=":")

            ax_amp_v_per.axhline(
                amp_ir_50,
                color="r",
                linestyle="-",
                lw=0.75,
                alpha=0.5,
                xmax=1.45,
                clip_on=False,
                zorder=100,
            )
            ax_amp_v_per.axhline(
                amp_nir_50,
                color="k",
                linestyle="--",
                lw=0.75,
                alpha=0.5,
                xmax=1.45,
                clip_on=False,
            )
            ax_amp_v_per.axvline(
                per_ir_50,
                color="r",
                linestyle="-",
                lw=0.75,
                alpha=0.5,
                ymax=1.45,
                clip_on=False,
            )
            ax_amp_v_per.axvline(
                per_nir_50,
                color="k",
                linestyle="--",
                lw=0.75,
                alpha=0.5,
                ymax=1.45,
                clip_on=False,
            )

            # fig.savefig()
            # return fig


# if __name__ == "__main__":

#     make_the_figure()
