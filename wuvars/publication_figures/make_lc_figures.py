"""
This script will do three things:

1. generate detrended, period-folded eight-panel lightcurves for all periodic variables
   (for FigureSet)
2. generate simple five-panel lightcurves for all non-periodic variables
3. generate simple five-panel lightcurves for all non-variables (for the appendix)

"""

import os
import pdb

import astropy.table
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.detrending import poly_detrend
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.data.preferred_photometry import (photometry_wserv7,
                                              photometry_wserv8)
from wuvars.plotting.lightcurve import (eightpanel_lc, eightpanel_lc_v2,
                                        ic348_eightpanel_lc,
                                        ic348_eightpanel_lc_v2, ic_date_offset,
                                        ngc1333_eightpanel_lc,
                                        ngc1333_eightpanel_lc_v2,
                                        ngc_date_offset)
from wuvars.publication_figures.make_figures import figure_export_path

ngc_match = match_ngc()
ic_match = match_ic()


region_keys = ["ic", "ngc"]
region_shorthand = {"ic": "ic348", "ngc": "ngc1333"}
wserv_dict = {"ngc": 7, "ic": 8}
fullname_dict = {"ngc": "NGC 1333", "ic": "IC 348"}
match_dict = {"ngc": ngc_match, "ic": ic_match}
index_dict = {"ic": "L16_T1_index", "ngc": "L16_T2_index"}
spread_dict = {"ngc": spreadsheet.load_wserv_v2(7), "ic": spreadsheet.load_wserv_v2(8)}
q_dict = {"ngc": quality_classes.load_q(7), "ic": quality_classes.load_q(8)}
data_dict = {
    "ngc": photometry_wserv7,
    "ic": photometry_wserv8,
}

detrended_data_path = "/Users/tsrice/Documents/Variability_Project_2020/Results/detrended_periodic_photometry"
detrended_data_dict = {
    "ngc": astropy.table.Table.read(
        os.path.join(
            detrended_data_path, "periodic_variables_detrended_photometry_ngc1333.h5"
        )
    ).group_by("SOURCEID"),
    "ic": astropy.table.Table.read(
        os.path.join(
            detrended_data_path, "periodic_variables_detrended_photometry_ic348.h5"
        )
    ).group_by("SOURCEID"),
}


lc_fn_dict = {"ngc": ngc1333_eightpanel_lc, "ic": ic348_eightpanel_lc}
lc_fn_dict_v2 = {"ngc": ngc1333_eightpanel_lc_v2, "ic": ic348_eightpanel_lc_v2}

# load up that ecsv
variability_tables = {
    "ngc": astropy.table.Table.read(
        os.path.join(figure_export_path, "table2_variability_ngc_ecsv.ecsv")
    ),
    "ic": astropy.table.Table.read(
        os.path.join(figure_export_path, "table2_variability_ic_ecsv.ecsv")
    ),
}
# cmap_dict = {"ngc": "jet", "ic": "jet_r"}

cmap_dict = {"ngc": "cubehelix_r", "ic": "cubehelix_r"}


date_offset_dict = {"ngc": ngc_date_offset, "ic": ic_date_offset}
date_breaks_dict = {"ngc": (0, np.inf), "ic": (60, 250)}


# make a periodic figure
def make_periodic_lc_figure(region_key, row, plotting=False, **kwargs):

    print(f"We're making a periodic figure in {region_key}")
    sid = row["SourceID"]
    name = row["Name"]
    period = row["Period"]
    detrend = row["PeriodDetrendMethod"]
    bestband = row["PeriodBand"]

    print(f"{sid = !s}, {name = !s}, {period=:.3f}, {bestband = !s}, {detrend = !s}")

    # okay! so in principle we have everything we need to make the figure...

    if plotting:

        dg = data_dict[region_key]

        if detrend == "vanilla":
            poly_order = 0
        else:
            poly_order = int(detrend[-1])

        dat_sid = dg.groups[dg.groups.keys["SOURCEID"] == sid]
        detrended, fit, fit_params_dict = poly_detrend(
            dat_sid,
            date_offset_dict[region_key],
            data_start=date_breaks_dict[region_key][0],
            data_end=date_breaks_dict[region_key][1],
            poly_order=poly_order,
            normalize=True,
            extrapolate_detrend=False,
            return_fit_params=True,
        )

        # This is where we would do some detrend stuff.

        with mpl.rc_context(
            {
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.top": True,
                "ytick.right": True,
            }
        ):
            fig = lc_fn_dict[region_key](
                dg,
                sid,
                period,
                detrended_dg=detrended,
                fit_dg=fit,
                fit_params_dict=fit_params_dict,
                **kwargs,
            )

            # EXTREMELY hacky workaround for brokenaxes
            for ax_band in [fig.ax_j, fig.ax_h, fig.ax_k]:
                for ax in ax_band.axs:

                    ax.xaxis.set_ticks_position("both")
                    # ax.xaxis.set_ticklabels([])

                    if len(ax_band.axs) == 1:
                        # ax.yaxis.set_ticks_position('both')
                        ax.tick_params(axis="y", left=True, right=True)
                        ax.tick_params(axis="y", labelright=False)

                    elif ax is ax_band.axs[-1]:
                        ax.yaxis.set_ticks_position("right")
                        ax.yaxis.set_ticklabels([])

        """
            fig_lc = ic348_simple_lc_scatter_brokenaxes(dat, sid, cmap='jet_r')    
            fig_lc = ngc1333_simple_lc_scatter_brokenaxes(dat, sid, cmap='jet')
        """

        plt.show()

    else:
        return None

    return fig


# make a periodic figure
def make_periodic_lc_figure_v2(region_key, row, plotting=False, **kwargs):

    sid = row["SourceID"]
    name = row["Name"]
    period = row["Period"]
    detrend = row["PeriodDetrendMethod"]
    bestband = row["PeriodBand"]

    print(f"{sid = !s}, {name = !s}, {period=:.3f}, {bestband = !s}, {detrend = !s}")

    # okay! so in principle we have everything we need to make the figure...

    if plotting:

        dg = detrended_data_dict[region_key]

        if detrend == "vanilla":
            poly_order = 0
        else:
            poly_order = int(detrend[-1])

        dat_sid = dg.groups[dg.groups.keys["SOURCEID"] == sid]

        # This is where we would do some detrend stuff.

        with mpl.rc_context(
            {
                "xtick.direction": "in",
                "ytick.direction": "in",
                "xtick.top": True,
                "ytick.right": True,
                "mathtext.fontset": "stix",       # use STIX fonts (Times-compatible)
                # "font.family": "serif",
                # "font.serif": ["Times New Roman"],  # or ["Times"]
            }
        ):
            fig = lc_fn_dict_v2[region_key](
                dg,
                sid,
                **kwargs,
            )

            # fig.suptitle(f"{sid = !s}, {name = !s}, {period=:.3f}, {bestband = !s}, {detrend = !s}")

            fig.ax_j.set_title(f"Short name: {name}, WFCAM ID: {sid}")
            fig.ax_j_phase.set_title(f"Best band: {bestband}, Detrend method: {detrend}")

            # EXTREMELY hacky workaround for brokenaxes
            for ax_band in [fig.ax_j, fig.ax_h, fig.ax_k]:
                for ax in ax_band.axs:

                    ax.xaxis.set_ticks_position("both")
                    # ax.xaxis.set_ticklabels([])

                    if len(ax_band.axs) == 1:
                        # ax.yaxis.set_ticks_position('both')
                        ax.tick_params(axis="y", left=True, right=True)
                        ax.tick_params(axis="y", labelright=False)

                    elif ax is ax_band.axs[-1]:
                        ax.yaxis.set_ticks_position("right")
                        ax.yaxis.set_ticklabels([])

        """
            fig_lc = ic348_simple_lc_scatter_brokenaxes(dat, sid, cmap='jet_r')    
            fig_lc = ngc1333_simple_lc_scatter_brokenaxes(dat, sid, cmap='jet')
        """


    else:
        return None

    return fig


def make_periodic_lc_figures(sample_only=True):
    """
    Makes the periodic light curve figures.

    By default, only generates the first one and/or a chosen sample
    (to be implemented).

    """

    # for each region...
    for region_key in region_keys:

        print("")
        print(f"Making periodic lc figures in {fullname_dict[region_key]}")
        print("")

        figureset_path_per = os.path.join(
            figure_export_path, f"Figure7Set_{region_key}"
        )
        os.makedirs(figureset_path_per, exist_ok=True)

        # grab the table.
        variability_table = variability_tables[region_key]

        periodic = variability_table["Periodic"] == "Y"
        variability_table_per = variability_table[periodic]

        # loop over the table
        for i, row in enumerate(variability_table_per):

            # if i<10:
            #     continue

            saving = True

            print(f"{i=}")


            # if i <= 30:
            #     continue
            # if i > 3:
            #     break

            # print the name of the source, its SOURCEID,
            # its period, and the detrending method.
            # Plot it. First in the stupidest possible way,
            #  but eventually with the cool fancy aspect.
            fig = make_periodic_lc_figure_v2(
                region_key, row, plotting=True, cmap=cmap_dict[region_key]
            )

            filename = os.path.join(figureset_path_per, f"{row['Name']}_lc")
            filename_pdf = filename + ".pdf"
            filename_png = filename + ".png"
            print(f"Saving fig to... {filename_pdf}")
            print(f"Saving fig to... {filename_png}")

            if saving:
                fig.savefig(filename_pdf, bbox_inches="tight")
                fig.savefig(filename_png, bbox_inches="tight")
                plt.close(fig)

    pass


def make_nonperiodic_variable_lc_figures(sample_only=True):

    pass


def make_nonvariable_lc_figures(sample_only=True):

    pass


if __name__ == "__main__":

    sample_only = True

    make_periodic_lc_figures(sample_only)
    make_nonperiodic_variable_lc_figures(sample_only)
    make_nonvariable_lc_figures(sample_only)
