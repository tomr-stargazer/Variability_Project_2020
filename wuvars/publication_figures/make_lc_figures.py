"""
This script will do three things:

1. generate detrended, period-folded eight-panel lightcurves for all periodic variables
   (for FigureSet)
2. generate simple five-panel lightcurves for all non-periodic variables
3. generate simple five-panel lightcurves for all non-variables (for the appendix)

"""

import os

import astropy.table
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.plotting.lightcurve import eightpanel_lc
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
    "ngc": photometry.group_wserv_v2(photometry.load_wserv_v2(7)),
    "ic": photometry.group_wserv_v2(photometry.load_wserv_v2(8))
}

# load up that ecsv
variability_tables = {
    "ngc": astropy.table.Table.read(
        os.path.join(figure_export_path, "table2_variability_ngc_ecsv.ecsv")
    ),
    "ic": astropy.table.Table.read(
        os.path.join(figure_export_path, "table2_variability_ic_ecsv.ecsv")
    ),
}


# make a periodic figure
def make_periodic_lc_figure(region_key, row, plotting=False):

    print(f"We're making a periodic figure in {region_key}")
    sid = row['SourceID']
    name = row['Name']
    period = row['Period']
    detrend = row['PeriodDetrendMethod']
    bestband = row['PeriodBand']

    print(f"{sid = !s}, {name = !s}, {period=:.3f}, {bestband = !s}, {detrend = !s}")

    # okay! so in principle we have everything we need to make the figure...

    if plotting:
        fig = eightpanel_lc()
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
        print(f"Doing thigns for {fullname_dict[region_key]}")
        print("") 


        figureset_path_per = os.path.join(
            figure_export_path, f"Figure7Set_{region_key}"
        )
        os.makedirs(figureset_path_per, exist_ok=True)

        # grab the table.
        variability_table = variability_tables[region_key]

        periodic = variability_table['Periodic'] == 'Y'
        variability_table_per = variability_table[periodic]
    
        # loop over the table
        for i, row in enumerate(variability_table_per):

            saving = False

            print(f"{i=}")

            if i>3: break


            # print the name of the source, its SOURCEID, 
            # its period, and the detrending method.
            # Plot it. First in the stupidest possible way, 
            #  but eventually with the cool fancy aspect.
            fig = make_periodic_lc_figure(region_key, row, plotting=True)

            filename = os.path.join(figureset_path_per, f"{row['Name']}_lc.pdf")
            print(f"Saving fig to... {filename}")

            if saving:
                fig.savefig(filename, bbox_inches='tight')

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
