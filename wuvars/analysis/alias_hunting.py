"""
Given a periodogram's peak frequency, how can we find its potential aliases?

Follows notation set out by vanderplas et al 2018
https://ui.adsabs.harvard.edu/abs/2018ApJS..236...16V/abstract

"""

import numpy as np


def find_aliases(f_peak, max_m=3, max_n=2):
    """
    Assume your highest peak is an alias of a true period. 
    What might the true period be?

    """

    return find_m_aliases(f_peak, max_m) + find_n_aliases(f_peak, max_n)


def find_m_aliases(f_peak, max_m=3):

    m_list = range(2, max_m + 1)

    f_aliases = []
    for m in m_list:
        for func in (np.multiply, np.divide):
            alias = func(f_peak, m)
            f_aliases.append(alias)

    return f_aliases


def find_n_aliases(f_peak, max_n=2):

    n_list = range(1, max_n + 1)

    f_aliases = []
    for n in n_list:
        for sign in (-1, +1):
            alias = np.abs(f_peak + sign * n)
            f_aliases.append(alias)

    return f_aliases


if __name__ == "__main__":

    np.testing.assert_equal(find_m_aliases(2.5, max_m=3), [5.0, 1.25, 7.5, 2.5 / 3])

    np.testing.assert_equal(find_n_aliases(2.5, max_n=2), [1.5, 3.5, 0.5, 4.5])
