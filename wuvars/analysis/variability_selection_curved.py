"""
Tools to select variable stars when the Stetson index declines vs. magnitude.

"""

import numpy as np
import pandas as pd
import pdb

def curve_Stetson(ds, magnitudes, thresholds, suffix, input_column='Stetson_JHK'):
    """
    Takes a summary spreadsheet `ds` and the output of  creates a new column with

    """

    k_mags = ds['median']['KAPERMAG3']
    stetson = ds['variability'][input_column]

    thresholds_interp = np.interp(k_mags, magnitudes, thresholds)

    subtracted = stetson - thresholds_interp

    ds['variability', input_column+suffix] = subtracted



def sv_jh(ds, Stetson_cutoff=3, suffix='', **kwargs):
    """
    Select 2-band variable stars, given a minimum JH Stetson index.

    """

    v_jh = ds["variability"]["Stetson_JH"+suffix] > Stetson_cutoff

    return v_jh


def sv_hk(ds, Stetson_cutoff=3, suffix='', **kwargs):
    """
    Select 2-band variable stars, given a minimum HK Stetson index.

    """

    v_hk = ds["variability"]["Stetson_HK"+suffix] > Stetson_cutoff

    return v_hk


def sv_jk(ds, Stetson_cutoff=3, suffix='', **kwargs):
    """
    Select 2-band variable stars, given a minimum JK Stetson index.

    """

    v_jk = ds["variability"]["Stetson_JK"+suffix] > Stetson_cutoff

    return v_jk


def sv_jhk(ds, Stetson_cutoff=5, suffix='', **kwargs):
    """
    Select 3-band variable stars, given a minimum JHK Stetson index.

    """

    v_jhk = ds["variability"]["Stetson_JHK"+suffix] > Stetson_cutoff

    return v_jhk


def sq2_variables(*args, **kwargs):
    """
    Select variable stars (double or triple-band) which meet Q2 quality standards.

    """
    ds = args[0]

    q2 = sq2(*args, **kwargs)

    v_jh = sv_jh(ds, **kwargs)
    v_hk = sv_hk(ds, **kwargs)
    v_jk = sv_jk(ds, **kwargs)

    v_jhk = sv_jhk(ds, **kwargs)

    q2_vars = q2 & (v_jhk | v_jk | v_hk | v_jh)

    return q2_vars


def sq1_variables(*args, **kwargs):
    """
    Select variable stars (double or triple-band) which meet Q1 quality standards.

    """
    ds = args[0]

    q1_j = sq1_j(*args, **kwargs)
    q1_h = sq1_h(*args, **kwargs)
    q1_k = sq1_k(*args, **kwargs)

    v_jh = sv_jh(ds, **kwargs)
    v_hk = sv_hk(ds, **kwargs)
    v_jk = sv_jk(ds, **kwargs)

    v_jhk = sv_jhk(ds, **kwargs)

    q1_vars = (
        (q1_j & q1_h & q1_k & v_jhk)
        | (q1_j & q1_h & v_jh)
        | (q1_j & q1_k & v_jk)
        | (q1_h & q1_k & v_hk)
    )

    return q1_vars


def sq0_variables(*args, **kwargs):
    """
    Select variable stars (single, double, or triple-band) which meet Q0 quality standards. 

    """

    ds = args[0]

    q0 = sq0(*args, **kwargs)

    v_j = sv_j(ds, **kwargs)
    v_h = sv_h(ds, **kwargs)
    v_k = sv_k(ds, **kwargs)

    v_jh = sv_jh(ds, **kwargs)
    v_hk = sv_hk(ds, **kwargs)
    v_jk = sv_jk(ds, **kwargs)

    v_jhk = sv_jhk(ds, **kwargs)

    q0_vars = q0 & (v_j | v_h | v_k | v_jh | v_jk | v_hk | v_jhk)

    return q0_vars


# This one was mainly used for exploratory purposes. See e.g., prototypes/Exploring K-only single band variables chisq.ipynb.
def sq1_k_variables_only(*args, **kwargs):
    ds = args[0]

    q1_j = sq1_j(*args, **kwargs)
    q1_h = sq1_h(*args, **kwargs)
    q1_k = sq1_k(*args, **kwargs)

    v_j = sv_j(ds, **kwargs)
    v_h = sv_h(ds, **kwargs)
    v_k = sv_k(ds, **kwargs)

    v_jh = sv_jh(ds, **kwargs)
    v_hk = sv_hk(ds, **kwargs)
    v_jk = sv_jk(ds, **kwargs)

    v_jhk = sv_jhk(ds, **kwargs)

    q1_k_only_vars = (
        (q1_k)
        & (v_k)
        & ~(
            (q1_j & v_j)
            | (q1_h & v_h)
            | (q1_j & q1_h & v_jh)
            | (q1_j & q1_k & v_jk)
            | (q1_h & q1_k & v_hk)
            | (q1_j & q1_h & q1_k & v_jhk)
        )
    )

    return q1_k_only_vars