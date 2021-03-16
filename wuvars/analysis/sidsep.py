"""
A short function that turns a 14-digit sourceid into a readable set of number clusters.

I don't plan to need much else here.

"""

import textwrap


def sidsep(sid):
    """
    Return a stringified sourceid with spaces.

    Turns
    44027709751316
    into
    4402 7709 7513 16
    .

    """

    s = str(sid)

    sep_sid = " ".join(textwrap.wrap(s, 4))

    return sep_sid
