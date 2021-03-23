"""
For now, this is a basic reduction placeholder script.

I want it to be a command-line tool that's like, 
>> python basic_reduction_script.py -i input_filename -o output_filename

and it just DOES IT.

Alternate usage:
(inside IPython/Jupyter)
%run basic_reduction_script.py # some arguments? not sure what's best

The following independent steps will occur in this:

- correct error bars
- clean the data by flagging and removing "bad" nights
- - first, find the constant stars: 
- - - run "naive" variability statistics over the whole dataset.
- - - select the stars that fall below an outlier-clipped median of naive-Stetson values, and a host of other criteria.
- - use those stars to "grade" each night
- - clean the data by removing nights below specific grades
- scrub data for other problematic datapoints
- - among stars with more than 90% good (unflagged) data, scrub out all flagged datapoints. (do this independently for each band)
- - remove all datapoints from the dataset that have photometric uncertainty greater than 0.5 mag.

Note: this is all following, in lock-step, the "detailed steps for processing of WFCAM data: Orion"

Some variations based on how the NGC 1333 data were reduced may arise. Maybe some command-line options will show up.

"""

import os
import numpy as np
import astropy.table

if __name__ == "__main__":

    print("Hello, world!")

    # just doing some prototyping here.

    data_location = os.path.expanduser(
        "~/Documents/Variability_Project_2020/wuvars/Data/Raw_Downloads/"
    )

    filename = "wserv8.fits.gz"
    filename_without_prefix = filename[: -len(".fits.gz")]
    print(filename_without_prefix)

    data = astropy.table.Table.read(data_location + filename)

    # some constants to manipulate the data
    s = 0.021
    c = 1.082

    data["JAPERMAG3ERR"] = np.sqrt(c * data["JAPERMAG3ERR"] ** 2 + s ** 2)
    data["HAPERMAG3ERR"] = np.sqrt(c * data["HAPERMAG3ERR"] ** 2 + s ** 2)
    data["KAPERMAG3ERR"] = np.sqrt(c * data["KAPERMAG3ERR"] ** 2 + s ** 2)

    data["JMHPNTERR"] = np.sqrt(c * data["JMHPNTERR"] ** 2 + s ** 2)
    data["HMKPNTERR"] = np.sqrt(c * data["HMKPNTERR"] ** 2 + s ** 2)

    target_location = os.path.expanduser(
        "~/Documents/Variability_Project_2020/wuvars/Data/reduction_artifacts/"
    )

    data.write(
        target_location + filename_without_prefix + "_revised_errorbars.fits",
        overwrite=True,
    )
    print(
        "Wrote new file to "
        + target_location
        + filename_without_prefix
        + "_revised_errorbars.fits"
    )

    # that should just about do it?
