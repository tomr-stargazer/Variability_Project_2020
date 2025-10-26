"""
UPDATE 26 October 2025: I think this code is deprecated and I am using make_lc_figures instead.

Okay. The whole point of this script is to generate some figures.

"""

import os

import astropy.table
import matplotlib.pyplot as plt
import numpy as np
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.data import photometry, quality_classes, spreadsheet
from wuvars.publication_figures.make_figures import figure_export_path

print("figureset_variables is deprecated.")

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

variability_tables = {
    "ngc": astropy.table.Table.read(
        os.path.join(figure_export_path, "table2_variability_ngc_ecsv.ecsv")
    ),
    "ic": astropy.table.Table.read(
        os.path.join(figure_export_path, "table2_variability_ic_ecsv.ecsv")
    ),
}


def make_nonperiodic_lc_figure(region_key, row):

    print("Nonperiodic figure NOT IMPLEMENTED YET")

    pass

# make a periodic figure
def make_periodic_lc_figure(region_key, row):

    print(f"We're making a periodic figure in {region_key}")
    sid = row['SourceID']
    name = row['Name']
    period = row['Period']
    detrend = row['PeriodDetrendMethod']
    bestband = row['PeriodBand']

    print(f"{sid = !s}, {name = !s}, {period=:.3f}, {bestband = !s}, {detrend = !s}")

    # okay! so in principle we have everything we need to make the figure...

    fig = eight

    return fig


def make_lc_figure_sets(ngc_posterboy=0, ic_posterboy=0):

    # load the IC 348 variability table

    # load the NGC 1333 variability table

    for region_key in region_keys:

        print("")
        print(f"Doing thigns for {fullname_dict[region_key]}")
        print("")

        figureset_path_per = os.path.join(
            figure_export_path, f"Figure7Set_{region_key}"
        )
        figureset_path_nonper = os.path.join(
            figure_export_path, f"Figure8Set_{region_key}"
        )

        os.makedirs(figureset_path_per, exist_ok=True)

        variability_table = variability_tables[region_key]

        for i, row in enumerate(variability_table):

            saving = False

            print(f"{i=}")

            if i>10: break
            
            if row['Periodic'] == 'Y':

                fig = make_periodic_lc_figure(region_key, row)

                filename = os.path.join(figureset_path_per, f"{row['Name']}_lc.pdf")
                print(f"Saving fig to... {filename}")
                if saving:
                    fig.savefig(filename, bbox_inches='tight')

            elif row['Periodic'] == 'N':

                fig = make_nonperiodic_lc_figure(region_key, row)

                filename = os.path.join(figureset_path_nonper, f"{row['Name']}_lc.pdf")
                print(f"Saving fig to... {filename}")
                if saving:
                    fig.savefig(filename, bbox_inches='tight')

            else:

                raise ValueError("Why is this row neither periodic nor nonperiodic?")

if __name__ == "__main__":

    make_lc_figure_sets()
