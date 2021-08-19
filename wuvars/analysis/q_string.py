"""
A function to make a little string which says what Q class an object is.

"""


def q_string(sid, spread, qualityset):

    q = qualityset

    if sid in spread[q.q2].index:
        return "2"

    elif sid in spread[q.q1_j | q.q1_h | q.q1_k].index:
        return_string = "1"

        if sid in spread[q.q1_j].index:
            return_string += "J"
        if sid in spread[q.q1_h].index:
            return_string += "H"
        if sid in spread[q.q1_k].index:
            return_string += "K"

        return return_string

    elif sid in spread[q.q0].index:
        return "0"

    else:
        return "-1"
