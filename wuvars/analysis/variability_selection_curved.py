"""
Tools to select variable stars when the Stetson index declines vs. magnitude.

"""

import numpy as np
import pandas as pd

def curve_Stetson(ds, magnitudes, thresholds, name=''):
    """
    Takes a summary spreadsheet `ds` and the output of  creates a new column with

    """

