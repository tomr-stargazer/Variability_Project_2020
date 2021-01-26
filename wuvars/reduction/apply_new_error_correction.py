"""
A function to revise errorbars.

There are two conceivable approaches:
- a single function that takes a filename, does the thing, writes a file, returns nothing
- two functions:
-   a function that takes just an astropy table and either alters it in-place or returns a revised one
-   a function that wraps the above fn with file read and write

The second approach is more testable. But this is SO simple.

"""

def revise_errorbars(improperly_error_corrected_dataset, s, c):

    # undo the old error correction terms

    # apply the new error correction terms

    # return the astropy table