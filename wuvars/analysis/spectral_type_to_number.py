"""
This is a quick little utility that turns a stringy spectral type into a number.

And vice versa.

It's really only intended for M and L types.

"""

import numpy as np


def get_num_from_SpT(SpT):
    """
    Turns a Spectral Type into a number.

    M0   = 0,
    M5.5 = 5.5,
    L1   = 11.

    """

    letter_part = SpT[0]

    letter_values = {
        "M": 0,
        "L": 10,
        "T": 20,
        "Y": 30,
        "K": -10,
        "G": -20,
        "F": -30,
        "A": -40,
        "B": -50,
        "O": -60,
    }  # newly extended to cover all the spectral types I think I might ever encounter

    number_part = float(SpT[1:])

    num = letter_values[letter_part] + number_part

    return num


def get_SpT_from_num(num):
    """
    Turns a number into a Spectral Type.

    M0   = 0,
    M5.5 = 5.5,
    L1   = 11.

    """

    letter_values = {
        "M": 0,
        "L": 10,
        "T": 20,
        "Y": 30,
        "K": -10,
        "G": -20,
        "F": -30,
        "A": -40,
        "B": -50,
        "O": -60,
    }  # newly extended to cover all the spectral types I think I might ever encounter

    letter_values_inv = {v: k for k, v in letter_values.items()}

    ones_place = num % 10
    tens_place = num - ones_place

    if isinstance(num, int):
        return f"{letter_values_inv[tens_place]}{ones_place}"
    elif isinstance(num, float):
        return f"{letter_values_inv[tens_place]}{ones_place:.1f}"


if __name__ == "__main__":

    assert get_num_from_SpT("M0") == 0
    assert get_num_from_SpT("M5.5") == 5.5
    assert get_num_from_SpT("L1") == 11
    assert get_num_from_SpT("K9") == -1
    assert get_num_from_SpT("G2") == -18

    assert get_SpT_from_num(0) == "M0"
    assert get_SpT_from_num(5.5) == "M5.5"
    assert get_SpT_from_num(11) == "L1"
    assert get_SpT_from_num(-1) == "K9"
    assert get_SpT_from_num(-18) == "G2"
