"""
Given a periodogram, how

"""

import numpy as np


# def forward_aliases(f_peak, max_m=3, max_n=2):
#     """
#     Assume your highest peak is the true frequency. 
#     What might its common aliases be?

#     """

#     pass

def find_aliases(f_peak, max_m=3, max_n=2):
    """
    Assume your highest peak is an alias of a true period. 
    What might the true period be?

    """

    return find_m_aliases(f_peak, max_m) + find_n_aliases(f_peak, max_n)

    # m_list = range(2, max_m + 1)
    # n_list = range(1, max_n + 1)

    # f_aliases = []

    # # for all combinations of m, n....
    # for m in m_list:
    #     print("m=" + str(m))
    #     for func in (np.multiply, np.divide):
    #         alias = func(f_peak, m)
    #         f_aliases.append(alias)

    # for n in n_list:
    #     print("n=" + str(n))
    #     for sign in (-1, +1):
    #         alias = np.abs(f_peak + sign * n)
    #         f_aliases.append(alias)

    # return f_aliases


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
