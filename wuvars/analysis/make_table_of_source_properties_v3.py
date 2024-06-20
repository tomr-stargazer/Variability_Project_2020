"""
This script will make a table of source properties - something that we may, later,
fill in (manually) with `variable? Y/N`, `periodic? Y/N`, etc.

Currently prototyping.

Columns I think I want, currently:
- Target name / ID
- RA/Dec
- median JHK
- Category (approved/statistical/color, or A/S/C)
- Variable? Y/N
- Type of variability (qualitative)
- Periodic? Y/N
- - Period
- - Amplitude (rms)
- delta J (90-10)
- delta H (90-10)
- delta K (90-10)

Something like that. Export it to a CSV or something.

So, for now, we want canonical source names etc.

"""

import numpy as np
import astropy.table

from wuvars.data import spreadsheet, photometry, quality_classes
from wuvars.analysis.spectral_type_to_number import get_num_from_SpT, get_SpT_from_num
from wuvars.analysis.bd_matching_v2 import match_onc, match_ngc, match_ic
from wuvars.analysis.q_string import q_string

# ngc_match = match_ngc()
# ic_match = match_ic()


def canonical_sourcenames(match, spread, qset, region: str):

    # build a string with the following:
    # which region
    # its running number among 'approved' sources (00-len-1)
    # a/s/c
    # teff or spectral class
    # q level
    #

    names = []

    for i, sid in enumerate(match.approved["SOURCEID"]):

        if sid in match.color["SOURCEID"]:
            asc_char = "C"
        elif sid in match.statistical["SOURCEID"]:
            asc_char = "S"
        else:
            asc_char = "A"

        if region == "ONC":
            status_str = f"{match.approved['Teff'][i]:4d}K"
        else:
            status_str = f"{get_SpT_from_num(match.approved['SpT'][i])}"

        qs = "Q" + q_string(sid, spread, qset)
        name = f"{region}_{i:03d}{asc_char}_{status_str}_{qs}"

        #         print(name)
        names.append(name)

    return names


def source_properties(matchstruct, spread, qset, region: str):
    """
    Columns I think I want, currently:
    - Target name / ID
    - RA/Dec
    - median JHK
    - Category (approved/statistical/color, or A/S/C)
    - Teff / SpT
    - Variable? Y/N
    - Type of variability (qualitative)
    - Periodic? Y/N
    - - Period
    - - Amplitude (rms)
    - delta J (90-10)
    - delta H (90-10)
    - delta K (90-10)    

    """

    table = astropy.table.Table()

    ma = matchstruct.approved

    names = canonical_sourcenames(matchstruct, spread, qset, region)
    asc_chars = []

    for i, sid in enumerate(matchstruct.approved["SOURCEID"]):

        if sid in matchstruct.color["SOURCEID"]:
            asc_chars.append("C")
        elif sid in matchstruct.statistical["SOURCEID"]:
            asc_chars.append("S")
        else:
            asc_chars.append("A")

    table["index"] = np.arange(len(ma["SOURCEID"]))
    table["SOURCEID"] = ma["SOURCEID"]
    table["shortname"] = names
    table["A/S/C"] = asc_chars
    table["RA_deg"] = np.degrees(ma["median_RA"])
    table["DE_deg"] = np.degrees(ma["median_DEC"])

    table["J_mag"] = ma["median_JAPERMAG3"]
    table["H_mag"] = ma["median_HAPERMAG3"]
    table["K_mag"] = ma["median_KAPERMAG3"]

    table["Teff"] = ma["Teff"]
    table["SpT"] = ma["SpT"]

    table["Stetson_JHK"] = ma["var_Stetson_JHK"]
    table["range_J"] = ma["range_JAPERMAG3"]
    table["range_H"] = ma["range_HAPERMAG3"]
    table["range_K"] = ma["range_KAPERMAG3"]

    table["period_J"] = np.zeros(len(table))
    table["per_amp_J"] = np.zeros(len(table))
    table["per_fap_J"] = np.zeros(len(table))

    table["period_H"] = np.zeros(len(table))
    table["per_amp_H"] = np.zeros(len(table))
    table["per_fap_H"] = np.zeros(len(table))

    table["period_K"] = np.zeros(len(table))
    table["per_amp_K"] = np.zeros(len(table))
    table["per_fap_K"] = np.zeros(len(table))

    return table


if __name__ == "__main__":

    onc_match = match_onc()
    ngc_match = match_ngc()
    ic_match = match_ic()

    onc_q = quality_classes.load_q(5)
    onc_spread = spreadsheet.load_wserv_v2(5)
    ngc_q = quality_classes.load_q(7)
    ngc_spread = spreadsheet.load_wserv_v2(7)
    ic_q = quality_classes.load_q(8)
    ic_spread = spreadsheet.load_wserv_v2(8)

    onc_table = source_properties(onc_match, onc_spread, onc_q, "onc")
    onc_table.write("ONC_source_properties.csv")

    ngc_table = source_properties(ngc_match, ngc_spread, ngc_q, "ngc")
    ngc_table.write("NGC_source_properties.csv")

    ic_table = source_properties(ic_match, ic_spread, ic_q, "ic")
    ic_table.write("IC_source_properties.csv")
