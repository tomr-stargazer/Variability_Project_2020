"""
This is a place to import data from. Here's where we can grab the photometry, 
as data structures.

Desired use:

from wuvars.data.photometry import v1

"""

import os
from datetime import datetime
from recordclass import recordclass
from astropy.table import Table

from wuvars.analysis.variability_selection import data_nuller

wserv_ids = [1, 5, 7, 8, 11]

wserv_ids_v2 = [5, 7, 8, 11]

Dataset = recordclass(
    "Dataset", (f"wserv{n}" for n in wserv_ids), defaults=(None,) * len(wserv_ids)
)

Dataset_v2 = recordclass(
    "Dataset", (f"wserv{n}" for n in wserv_ids_v2), defaults=(None,) * len(wserv_ids_v2)
)


def load_v1():

    v1 = Dataset()

    data_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"

    for i, wserv in enumerate(wserv_ids):
        data_path = os.path.join(
            data_root,
            f"wserv{str(wserv)}",
            f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_new_error_corrected.fits",
        )

        print(f"Loading WSERV{wserv} photometry data... ", end="", flush=True)
        startTime = datetime.now()

        dat = Table.read(data_path)
        v1[i] = dat

        print(
            f"DONE (elapsed time: {(datetime.now() - startTime).total_seconds():.2f}s)"
        )

    return v1


def load_v1_grouped(v1=None):

    if v1 is None:
        v1 = load_v1()

    v1_grouped = Dataset()

    for i, wserv in enumerate(wserv_ids):

        print(f"Grouping WSERV{wserv} photometry data... ", end="", flush=True)
        startTime = datetime.now()

        dat = v1[i]
        df = dat.to_pandas()
        data_nuller(df)
        dat_again = Table.from_pandas(df)
        dat_by_source = dat_again.group_by("SOURCEID")
        v1_grouped[i] = dat_by_source

        print(
            f"DONE (elapsed time: {(datetime.now() - startTime).total_seconds():.2f}s)"
        )

    return v1_grouped


# load and group for each


def load_wserv_v2(wserv):
    """
    Loads a photometry data given a WSERV id.

    """
    data_root = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"
    data_path = os.path.join(
        data_root,
        f"wserv{str(wserv)}",
        f"WSERV{str(wserv)}_graded_clipped0.95_scrubbed0.1_dusted0.5_new_error_corrected.fits",
    )

    print(f"Loading WSERV{wserv} photometry data... ", end="", flush=True)
    startTime = datetime.now()

    dat = Table.read(data_path)
    print(f"DONE (elapsed time: {(datetime.now() - startTime).total_seconds():.2f}s)")
    return dat


def group_wserv_v2(dat):

    df = dat.to_pandas()
    data_nuller(df)
    dat_again = Table.from_pandas(df)
    dat_by_source = dat_again.group_by("SOURCEID")

    return dat_by_source


def load_v2():
    v2 = Dataset_v2()
    for i, wserv in enumerate(wserv_ids_v2):
        v2[i] = load_wserv_v2(wserv)

    return v2


def load_v2_grouped(v2=None):

    if v2 is None:
        v2 = load_v2()

    v2_grouped = Dataset_v2()

    for i, wserv in enumerate(wserv_ids_v2):

        print(f"Grouping WSERV{wserv} photometry data... ", end="", flush=True)
        startTime = datetime.now()

        dat = v2[i]
        v2_grouped[i] = group_wserv_v2(dat)

        print(
            f"DONE (elapsed time: {(datetime.now() - startTime).total_seconds():.2f}s)"
        )

    return v2_grouped
