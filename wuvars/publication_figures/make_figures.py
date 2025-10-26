"""
This is a script designed for one thing: make all figures.

The goal is for it to be up to date and ease my process. 
Ideally: every publication figure (including appendices) is made here.

"""


figure_export_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/figs_tables_for_overleaf"
figureset_export_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/figure_sets"

if __name__ == "__main__":

    # Turn on when we have the full figure set making function
    if False:
        full_set = (input("Make full figure set? [y/N]: ").strip().lower() == "y")

    # Figure 1:

    # Mass-SpT relationship

    from mass_SpT_figure import make_mass_SpT_figure

    make_mass_SpT_figure()


    # Figure 2:

    # The Map

    from map_figure_Perseus_regions_v4 import maps_v4_plus_WISEcolor

    maps_v4_plus_WISEcolor()

    # Figure 3 + 4:

    # It's doing double duty! There's a CMD also being generated.
    from HR_diagram_NGC_IC import make_fig_3a_3b_HR_diagram

    make_fig_3a_3b_HR_diagram()

    # Figure -4- 5:
    # Period injection results

    from period_injection_figure import make_period_injection_fig

    make_period_injection_fig()

    # Figure 5: 
    # Period comparison

    from wuvars.analysis.period_match import make_period_comparison_figure
    make_period_comparison_figure(verbose=False, saving=True)

    # Figure 6:
    # FigureSet: periodic light curves

    # Figure 7:
    # FigureSet: non-periodic, variable light curves
