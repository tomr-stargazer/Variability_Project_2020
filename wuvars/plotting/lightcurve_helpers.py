"""
Some auxiliary functions to help the lightcurves.

At the moment I am mostly thinking of a way to produce `xlims`
for use in `brokenaxes` axes, given arrays of date values.

I've also copied in the `orion_cmap` that I used for the 2015 paper. 
It's really useful and well-tuned to its use case.

"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap


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
    xlims = [(xlim_list[i], xlim_list[i + 1]) for i in range(0, len(xlim_list), 2)]

    return xlims


# fmt: off
orion_color_dic = {'red': ((0., 1, 0),
                           (0.4, 1, 1),
                           (0.9, 1, 1),
                           (1, 0.5, 0.5)),
                   'green': ((0., 1, 0),
                             (0.136, 1, 1),
                             (0.8, 1, 1),
                             (0.9, 0, 0),
                             (1, 0, 0)),
                   'blue': ((0., 1, 1),
                            (0.2, 0.5, 0.5),
                            (0.7, 0, 0),
                            (1, 0, 0))}
# fmt: on

orion_cmap = LinearSegmentedColormap("orion_cmap", orion_color_dic)
