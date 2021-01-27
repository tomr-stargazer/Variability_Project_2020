"""
A function to revise errorbars.

There are two conceivable approaches:
- a single function that takes a filename, does the thing, writes a file, returns nothing
- two functions:
-   a function that takes just an astropy table and either alters it in-place or returns a revised one
-   a function that wraps the above fn with file read and write

The second approach is more testable. But this is SO simple.

"""

import numpy as np
from wuvars.analysis.rederive_error_correction import decorrect_error


def revise_errorbars(improperly_error_corrected_dataset, s, c):

    dat = improperly_error_corrected_dataset

    bands = ["J", "H", "K"]
    for b in bands:

        # undo the old error correction terms

        improper_err = dat[f"{b}APERMAG3ERR"]
        err = decorrect_error(improper_err)

        # apply the new error correction terms
        new_corrected_err = np.sqrt(c * err ** 2 + s ** 2)

        dat[f"{b}APERMAG3ERR"] = new_corrected_err

    return dat


if __name__ == "__main__":
    print("Testing")
