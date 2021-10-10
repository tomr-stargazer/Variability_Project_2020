"""
Tom is rewriting this around the following use-case.

We'll input
(a)  a list of literature catalogs of brown dwarfs, and
(b) our own spreadsheet,
and we'll expect as output
(a) a list (table?) of successful cross-matches, and 
(b) a list (table?) of failed cross-matches.


"""

import os
from collections import OrderedDict, namedtuple

import pdb

import numpy as np
import pandas as pd

import astropy.table
from astropy import units as u
from astropy.coordinates import SkyCoord


# Table_and_coord = namedtuple("Table_and_coord", ("table", "coords"))


class TableMatcher(object):
    """
    Basically, I'm implementing the idea I'd written up in my notes.

    """

    def __init__(self, spread, aux_tables):
        super(TableMatcher, self).__init__()
        self.spread = spread
        self.aux_tables = aux_tables
        self.match_table = None
        self.non_matches = None
        self.BD_spreadsheet = None

    def match(self):
        """
        We'll implement stuff here from my notes.

        """
        BD_match_sourceids_list = []
        non_matches_dict = {}

        for (name, aux_table) in self.aux_tables.items():

            idx, d2d, d3d = aux_table.coords.match_to_catalog_sky(self.spread.coords)

            max_sep = 1.0 * u.arcsec
            sep_constraint = d2d < max_sep

            aux_table.matches = self.spread.iloc[idx[sep_constraint]]
            aux_table.non_matches = aux_table[~sep_constraint]

            # print("Matches: \n", matches, "\n\n")
            print(f"Non-matches to {name}:", len(aux_table.non_matches))

            BD_match_sourceids_list.extend(aux_table.matches.index)
            non_matches_dict[name] = aux_table.non_matches

        BD_match_sourceids = np.unique(BD_match_sourceids_list)
        BD_spreadsheet = self.spread[np.in1d(self.spread.index, BD_match_sourceids)]
        BD_coords = self.spread.coords[np.in1d(self.spread.index, BD_match_sourceids)]

        match_table = astropy.table.Table()
        match_table["SOURCEID"] = BD_spreadsheet.index
        match_table["RA"] = np.degrees(BD_spreadsheet["median"]["RA"])
        match_table["DEC"] = np.degrees(BD_spreadsheet["median"]["DEC"])

        for (name, aux_table) in self.aux_tables.items():

            # note that the catalog matching here is intentionally REVERSE of previous
            idx, d2d, d3d = BD_coords.match_to_catalog_sky(aux_table.coords)

            max_sep = 1.0 * u.arcsec
            sep_constraint = d2d < max_sep

            # BD_matches = BD_coords[sep_constraint]
            # aux_matches = aux_table.coords[idx[sep_constraint]]

            mask_these_as_false = ~sep_constraint

            idx[~sep_constraint] = -99999
            idx_masked = np.ma.masked_array(idx, mask=mask_these_as_false)

            # pdb.set_trace()

            # print(np.max(d2d[sep_constraint].to(u.arcsec)))

            match_table[f"{name}_index"] = idx_masked

        # print("Match table:\n", match_table)

        self.match_table = match_table
        self.BD_spreadsheet = BD_spreadsheet
        self.non_matches = non_matches_dict

        return None


if __name__ == "__main__":

    # we'll prototype what goes here in a jupyter notebook.
    # I'm mostly thinking it'll be NGC 1333 as an example.

    aux_path = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/data/auxiliary_catalogs/NGC1333"

    filepath_table4 = os.path.join(aux_path, "Scholz_2012ApJ_744_6S_table4.fit")
    table4 = astropy.table.Table.read(filepath_table4)
    table4.coords = SkyCoord(ra=table4["RAJ2000"], dec=table4["DEJ2000"])
    # table4_tc = Table_and_coord(table4, coords_table4)

    filepath_table2 = os.path.join(aux_path, "Scholz_2012ApJ_744_6S_table2.fit")
    table2 = astropy.table.Table.read(filepath_table2)
    table2.coords = SkyCoord(ra=table2["RAJ2000"], dec=table2["DEJ2000"])
    # table2_tc = Table_and_coord(table2, coords_table2)

    # print(table2.coords)
    # print(table4.coords)

    filepath_table1 = os.path.join(aux_path, "Scholz_2012ApJ_756_24S_table1.fit")
    table1 = astropy.table.Table.read(filepath_table1)
    table1.coords = SkyCoord(ra=table1["RAJ2000"], dec=table1["DEJ2000"])
    # table1_tc = Table_and_coord(table1, coords_table1)

    # print(table1.coords)

    # this is a sufficient basis for me to now build my thing around.

    aux_tables = OrderedDict()
    aux_tables["Scholz_12a_table2"] = table2
    aux_tables["Scholz_12a_table4"] = table4
    aux_tables["Scholz_12b_table1"] = table1

    # now we grab our spreadsheet
    from wuvars.data import spreadsheet

    wserv = 7  # for NGC 1333
    spread = spreadsheet.load_wserv_v2(wserv)
    spread.coords = SkyCoord(
        ra=spread["median"]["RA"].values * u.rad,
        dec=spread["median"]["DEC"].values * u.rad,
    )

    # spread_tc = Table_and_coord(spread, spreadsheet_coordinates)

    if False:

        # I think that everything below here can go in a method or function which takes
        # the following as input:
        #    spread_tc:  Table_and_coord of (pandas.DataFrame, SkyCoord)
        #    aux_tables:  Table_and_coord of (astropy.Table, SkyCoord)
        #
        # and then "returns" (sets the variables in `self`) the following:
        #    BD_spread_tc:  Table_and_coord of (pandas.DataFrame, SkyCoord)
        #    match_table:

        # ngc1333 = TableMatcher(spread_tc, aux_tables)
        # ngc1333.match()
        # ngc1333.match_table
        # ngc1333.BD_spreadsheet
        # ngc1333.non_matches_dict

        BD_match_sourceids_list = []

        # non_matches_dict = {}

        for (name, aux_table) in aux_tables.items():

            idx, d2d, d3d = aux_table.coords.match_to_catalog_sky(spread.coords)

            max_sep = 1.0 * u.arcsec
            sep_constraint = d2d < max_sep

            aux_table.matches = spread.iloc[idx[sep_constraint]]
            aux_table.non_matches = aux_table[~sep_constraint]

            # print("Matches: \n", matches, "\n\n")
            # print("Non-matches: \n", non_matches, "\n\n")

            BD_match_sourceids_list.extend(aux_table.matches.index)
            # non_matches_dict[name] = non_matches

        BD_match_sourceids = np.unique(BD_match_sourceids_list)
        BD_spreadsheet = spread[np.in1d(spread.index, BD_match_sourceids)]
        BD_coords = spread.coords[np.in1d(spread.index, BD_match_sourceids)]

        match_table = astropy.table.Table()
        match_table["SOURCEID"] = BD_spreadsheet.index
        match_table["RA"] = np.degrees(BD_spreadsheet["median"]["RA"])
        match_table["DEC"] = np.degrees(BD_spreadsheet["median"]["DEC"])

        for (name, aux_table) in aux_tables.items():

            # note that the catalog matching here is intentionally REVERSE of previous
            idx, d2d, d3d = BD_coords.match_to_catalog_sky(aux_table.coords)

            max_sep = 1.0 * u.arcsec
            sep_constraint = d2d < max_sep

            # BD_matches = BD_coords[sep_constraint]
            # aux_matches = aux_tc.coords[idx[sep_constraint]]

            idx[~sep_constraint] = -99999

            print(np.max(d2d[sep_constraint].to(u.arcsec)))

            match_table[f"{name}_index"] = idx

        print("Match table:\n", match_table)
