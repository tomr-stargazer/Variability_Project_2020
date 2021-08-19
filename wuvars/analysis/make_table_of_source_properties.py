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
