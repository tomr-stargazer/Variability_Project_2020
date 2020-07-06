# This loads each dataset and does the simplest possible thing:
# it prints the number of rows.

import os
import glob

data_location = os.path.expanduser(
    "~/Documents/Variability_Project_2020/Data/Raw_Downloads/"
)

data_files = glob.glob(data_location + "*.fits.gz")

# zero-order sanity check: print the filenames with their paths.

print ("Here are the filenames: ")
for filename in data_files:
    print(filename)
print("")

# first, sanity check, get the actual filesizes in bytes
print("Here are their sizes, in bytes")
for filename in data_files:
    print(os.path.getsize(filename))

# # now, load each file, print its size (as a length), then close it
# for filename in data_files:

# 	# this `tab` variable is going to be overwritten in every iteration of the loop
# 	tab =
#     print(os.path.getsize(filename))
