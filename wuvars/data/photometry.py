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

Dataset = recordclass(
    "Dataset", (f"wserv{n}" for n in wserv_ids), defaults=(None,) * len(wserv_ids)
)


v1 = Dataset()

v1_grouped = Dataset()


data_root = (
    "/Users/tsrice/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"
)

photometry_data = []

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
    print("Grouping... ", end="", flush=True)

    df = dat.to_pandas()
    data_nuller(df)
    dat_again = Table.from_pandas(df)
    dat_by_source = dat_again.group_by("SOURCEID")
    v1_grouped[i] = dat_by_source

    print(f"DONE (elapsed time: {(datetime.now() - startTime).total_seconds():.2f}s)")
