# I want a summary plot showing the time coverage of each of the five surveys.
# It'll be simple, just versus MJD, nothing fancy.

# I'm doing imports following PEP8, as summarized here:
# https://stackoverflow.com/questions/22722976/import-order-coding-standard
import os
import glob
import time

import numpy as np
import astropy.table as tab
import matplotlib.pyplot as plt


data_location = os.path.expanduser(
    "~/Documents/Variability_Project_2020/Data/Raw_Downloads/"
)

data_files = glob.glob(data_location + "*.fits.gz")

print([x[len(data_location) :] for x in data_files])

fig = plt.figure(figsize=(10, 4))

# I want them to go in order of first observation. This is an ad-hoc ordering list.
top_to_bottom_order = (0, 3, 4, 1, 2)
plot_labels = ["WSERV1", "WSERV5", "WSERV7", "WSERV8", "WSERV11"]

for i, filename in enumerate([data_files[x] for x in top_to_bottom_order][::-1]):
    start_time = time.time()

    # load the data table
    table = tab.Table.read(filename, format="fits")

    # this is a completely arbitrary slice that is meant to make it MUCH faster to generate the desired list of timestamps
    # the slice along J magnitude should be completely orthogonal to dates, so in principle no information should be lost
    table_subset = table[(table["JAPERMAG3"] > 12) & (table["JAPERMAG3"] < 12.5)]

    # get all the dates in this table, throwing away duplicates
    dates = np.unique(table_subset["MEANMJDOBS"])
    print("number of dates: ", len(dates))

    plt.plot(dates, [i] * len(dates), "k.")
    elapsed_time = time.time() - start_time

    print(f"Elapsed time for {i}: {elapsed_time:.2f} s")

    # this table object is very memory-intensive, so we're going to delete it and hope it gets garbage collected.
    # This is not standard practice in Python programming; usually you don't need to delete variables when you're done.
    del table
    del table_subset

# get current axis (gca)
ax = plt.gca()
ax.set_xlabel("Modified Julian Date")
ax.set_yticks(np.arange(5))
ax.set_yticklabels(plot_labels[::-1], rotation="horizontal")

# some vertical lines marking january 1 of various years:
for i in range(11):  # goes from 0 to 10
    # MJD 54101 is January 1, 2007
    ax.axvline(54101 + 365.25 * i, color="k", alpha=0.5, lw=0.5)

plt.show()
fig.savefig("date_coverage_by_program.png")
