"""
This is a script designed for one thing: make all figures.

The goal is for it to be up to date and ease my process. 
Ideally: every publication figure (including appendices) is made here.

"""


figure_export_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/figs_tables_for_overleaf"


if __name__ == "__main__":

    # Figure 1:

    # Mass-SpT relationship

    from mass_SpT_figure import make_mass_SpT_figure

    make_mass_SpT_figure()


    # Figure 2:

    # The Map

    from map_figure_Perseus_regions_v4 import maps_v4_plus_WISEcolor

    maps_v4_plus_WISEcolor()

    # Figure 3:

    from HR_diagram_NGC_IC import make_fig_3a_3b_HR_diagram

    make_fig_3a_3b_HR_diagram()


    # (What is it?)
