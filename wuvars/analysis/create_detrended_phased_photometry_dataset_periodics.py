"""
So. This is a script that loops through all of the periodic variables,
does the detrending, and generates a new photometry dataset which
contains exactly all the columns of the original , PLUS

- period*
- trend_params*
- phase
- j_corr
- h_corr
- k_corr

* these columns are static re: source.

"""

import os
import pdb

import astropy.table
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import MaskedColumn
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.detrending import poly_detrend
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.data.preferred_photometry import (photometry_wserv7,
                                              photometry_wserv8)
from wuvars.plotting.lightcurve import (eightpanel_lc, ic348_eightpanel_lc,
                                        ic_date_offset, ngc1333_eightpanel_lc,
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
lc_fn_dict = {"ngc": ngc1333_eightpanel_lc, "ic": ic348_eightpanel_lc}

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

cmap_dict = {"ngc": "magma_r", "ic": "magma_r"}


date_offset_dict = {"ngc": ngc_date_offset, "ic": ic_date_offset}
date_breaks_dict = {"ngc": (0, np.inf), "ic": (60, 250)}

data_export_path = "/Users/tsrice/Documents/Variability_Project_2020/Results/detrended_periodic_photometry"


def make_detrended_photometry_dataset(region_key, verbose=False):
    """
    This takes a region_key and does all the work of updating the data,
    outputting a new "detrended" dataset.

    """

    dg = data_dict[region_key]

    variability_table = variability_tables[region_key]

    periodic = variability_table["Periodic"] == "Y"
    variability_table_per = variability_table[periodic]

    # OKAY maybe we should actually first restrict it to just the SOURCEIDs in the var
    return_dg = dg[np.in1d(dg["SOURCEID"], variability_table_per["SourceID"])].copy()
    blank_float_column = np.ones_like(return_dg["MEANMJDOBS"]) * np.nan
    blank_str_column = np.ones_like(return_dg["MEANMJDOBS"]).astype('str')

    new_column_names = [
        "period",
        "detrend",
        "bestband",
        "phase",
        "j_corr",
        "h_corr",
        "k_corr",
    ]

    poly_colnames = [f"{band}_poly_{n}" for band in ["J", "H", "K"] for n in range(5)]

    new_column_names.extend(poly_colnames)

    for col_name in new_column_names:
        if col_name in ['detrend', 'bestband']:
            return_dg.add_column(blank_str_column, name=col_name)
        elif 'corr' in col_name:
            band_up = col_name[0].upper()
            return_dg.add_column(MaskedColumn(blank_float_column, name=col_name, mask=return_dg[f"{band_up}APERMAG3"].mask))
        else:    
            return_dg.add_column(blank_float_column, name=col_name)

    # loop over the table
    for i, row in enumerate(variability_table_per):

        if verbose:
            print(f"{i=}")

        sid = row["SourceID"]
        name = row["Name"]
        period = row["Period"]
        detrend = row["PeriodDetrendMethod"]
        bestband = row["PeriodBand"]

        dat_sid = dg.groups[dg.groups.keys["SOURCEID"] == sid]

        if detrend == "vanilla":
            poly_order = 0
        else:
            poly_order = int(detrend[-1])

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

        # some data can be updated in bulk.
        bulk_update_selection = (return_dg["SOURCEID"] == sid)

        return_dg['period'][bulk_update_selection] = period
        return_dg['detrend'][bulk_update_selection] = detrend
        return_dg['bestband'][bulk_update_selection] = bestband
        for k in range(5):
            for band in ["J", "H", "K"]:
                try:
                    return_dg[f'{band}_poly_{k}'][bulk_update_selection] = fit_params_dict[band][k]

                except IndexError:
                    # usually because of a lower degree polynomial fit
                    pass
                except TypeError:
                    # usually because no data, so no polynomial fit?
                    pass
                

        # updating data parameters, one timestamp at a time
        for j, phot_row in enumerate(detrended):

            if verbose:
                print(f" {j=}")

            timestamp = phot_row["MEANMJDOBS"]

            update_row_selection = (return_dg["MEANMJDOBS"] == timestamp) & (
                return_dg["SOURCEID"] == sid
            )

            phase = (timestamp % period) / period

            return_dg['phase'][update_row_selection] = phase

            return_dg['j_corr'][update_row_selection] = phot_row['JAPERMAG3']
            return_dg['h_corr'][update_row_selection] = phot_row['HAPERMAG3']
            return_dg['k_corr'][update_row_selection] = phot_row['KAPERMAG3']


    # now apply masks
    # return_dg['j_corr'].mask = return_dg['JAPERMAG3'].mask
    # return_dg['h_corr'].mask = return_dg['HAPERMAG3'].mask
    # return_dg['k_corr'].mask = return_dg['KAPERMAG3'].mask

    if True:
        filename = os.path.join(
            data_export_path,
            f"periodic_variables_detrended_photometry_{region_shorthand[region_key]}.h5",
        )

        return_dg.write(filename, overwrite=True)
        if verbose:
            print(f"Wrote to {filename}")

        return return_dg


# make a periodic figure
def _make_periodic_lc_figure(region_key, row, plotting=False, **kwargs):

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


def _make_periodic_lc_figures(sample_only=True):
    """
    Makes the periodic light curve figures.

    By default, only generates the first one and/or a chosen sample
    (to be implemented).

    """

    # for each region...
    for region_key in region_keys:

        print("")
        print(f"Doing thigns for {fullname_dict[region_key]}")
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

            saving = False

            print(f"{i=}")

            if i > 3:
                break

            # print the name of the source, its SOURCEID,
            # its period, and the detrending method.
            # Plot it. First in the stupidest possible way,
            #  but eventually with the cool fancy aspect.
            fig = make_periodic_lc_figure(
                region_key, row, plotting=True, cmap=cmap_dict[region_key]
            )

            filename = os.path.join(figureset_path_per, f"{row['Name']}_lc.pdf")
            print(f"Saving fig to... {filename}")

            if saving:
                fig.savefig(filename, bbox_inches="tight")

    pass
