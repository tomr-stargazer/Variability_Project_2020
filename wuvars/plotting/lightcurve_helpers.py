"""
Some auxiliary functions to help the lightcurves.

At the moment I am mostly thinking of a way to produce `xlims`
for use in `brokenaxes` axes, given arrays of date values.

"""

import numpy as np


def produce_xlims(dates, breaks, pad=10):

    xlim_list = []
    xlim_list.append(np.floor(np.min(dates)) - pad)

    for _break in breaks:

        try:

            # left = np.floor(np.min(dates[dates < _break]))
            left = np.ceil(np.max(dates[dates < _break]) + pad)
            right = np.floor(np.min(dates[dates >= _break]) - pad)

            xlim_list.extend([left, right])

        except ValueError:
            pass

    # finally:
    xlim_list.append(np.ceil(np.max(dates)) + pad)

    # group xlim_list into pairs.
    xlims = [(xlim_list[i], xlim_list[i+1]) for i in range(0, len(xlim_list), 2)]

    return xlims
