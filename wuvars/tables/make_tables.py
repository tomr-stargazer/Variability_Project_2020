"""
This is a script designed for one thing: make all tables.

The goal is for it to be up to date and ease my process. 
Ideally: 

"""


table_export_path = (
    "/Users/tsrice/Documents/Variability_Project_2020/wuvars/figs_tables_for_overleaf"
)


if __name__ == "__main__":

    # Make the first table, in some script somewhere
    # We are going to (simultaneously in most cases)

    # Table(s) 1:

    from wuvars.tables.make_table_1_targets import make_table_1_targets

    make_table_1_targets(verbose=False)

    # Table(s) 2:

    from wuvars.tables.make_table_2_variability_properties import (
        make_table_2_variability_properties,
    )

    make_table_2_variability_properties(verbose=False)
