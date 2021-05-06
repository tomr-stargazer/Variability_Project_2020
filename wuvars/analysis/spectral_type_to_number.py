"""
This is a quick little utility that turns a stringy spectral type into a number.

And vice versa.

It's really only intended for M and L types.

"""


def get_num_from_SpT(SpT):
    """
    Turns a Spectral Type into a number.

    M0   = 0,
    M5.5 = 5.5,
    L1   = 11.

    """

    letter_part = SpT[0]

    letter_values = {"M": 0, "L": 10}

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

    if num < 10:
        if isinstance(num, int):
            return f"M{num}"
        elif isinstance(num, float):
            return f"M{num:.1f}"
    elif num >= 10:
        if isinstance(num, int):
            return f"L{num-10}"
        elif isinstance(num, float):
            return f"L{num-10:.1f}"


if __name__ == "__main__":

    assert get_num_from_SpT("M0") == 0
    assert get_num_from_SpT("M5.5") == 5.5
    assert get_num_from_SpT("L1") == 11

    assert get_SpT_from_num(0) == "M0"
    assert get_SpT_from_num(5.5) == "M5.5"
    assert get_SpT_from_num(11) == "L1"
