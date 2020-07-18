# This loads each dataset and does the simplest possible thing:
# it prints the number of rows.

import os
import glob
import time
import astropy.table as tab

data_location = os.path.expanduser(
    "~/Documents/Variability_Project_2020/wuvars/Data/Raw_Downloads/"
)

data_files = glob.glob(data_location + "*.fits.gz")

# zero-order sanity check: print the filenames with their paths.

print("\n  Here are the filenames: ")
for filename in data_files:
    print(filename)

data_files_without_prefix = [x[len(data_location) :] for x in data_files]

# first, sanity check, get the actual filesizes in bytes
print("\n  Here are their sizes, in MB:")
for filename, filename_sans_prefix in zip(data_files, data_files_without_prefix):
    print(f"{filename_sans_prefix}: {os.path.getsize(filename) / 1e6:.1f} MB")

# now, load each file, print its size (its length), then close it

print("\n  Load times and table sizes, in rows:")
for filename, filename_sans_prefix in zip(data_files, data_files_without_prefix):
    start_time = time.time()

    # this `table` variable is going to be overwritten in every iteration of the loop
    table = tab.Table.read(filename, format="fits")
    elapsed_time = time.time() - start_time

    print(f"{filename_sans_prefix} took {elapsed_time:.2f} s to load")
    print(f"{filename_sans_prefix} has {len(table)} rows ({len(table)/1e6:.1f}M rows)")
    print("")
