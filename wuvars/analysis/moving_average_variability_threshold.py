"""
Here I'm writing a script which will calculate a "moving average" variability threshold.

Here's the core of the idea:

The Stetson variability index has a different distribution at different magnitudes 
(its mean declines with fainter magnitudes). Therefore a "flat" cutoff is inappropriate;
it ignores demonstrably 

In order to estimate a sensible threshold...
- exclude the spatial region corresponding to where most of the YSOs are
- 

"""

import numpy as np


def moving_average_percentile(
    x_array,
    y_array,
    x_range,
    stepsize=0.1,
    percentile=95,
    window_size=0.75,
    ignore_y_above=2,
    verbose=False,
):
    """
    Takes in a dataset of (x, y) coordinates, 

    """

    # construct a grid
    x_grid = np.arange(*x_range, stepsize)
    result_grid = np.zeros_like(x_grid)

    for i, x in enumerate(x_grid):

        try:

            min_x = x - window_size
            max_x = x + window_size

            valid_ys = y_array[
                (x_array > min_x) & (x_array < max_x) & (y_array <= ignore_y_above)
            ]

            result = np.percentile(valid_ys, percentile)
            result_grid[i] = result
        except IndexError:
            result_grid[i] = np.nan
        if verbose:
            print(f"{x:.1f}, {result:.1f}")

    return x_grid, result_grid


def monotone_decline(xs_):

    xs = np.copy(xs_)

    for i in range(1, len(xs)):

        if xs[i] > xs[i - 1]:
            xs[i] = xs[i - 1]

    return xs


def moving_average_percentile_monotonic(*args, **kwargs):
    x_grid, result_grid = moving_average_percentile(*args, **kwargs)
    result_grid_monotonic = monotone_decline(result_grid)
    return x_grid, result_grid_monotonic
