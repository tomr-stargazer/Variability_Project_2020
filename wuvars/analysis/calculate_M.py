"""
I'm going to calculate M, the assymmetry metric identified by Cody et al. 2014.

I'm adopting the simplification described in Bredall+20 and Hillenbrand+22.

"""

import numpy as np

# Here's the main challenge:
# Segmenting the data into deciles!


def compute_M(dataset, sid):
    """ 
    Given a photometry dataset and a SID, I want you to compute M.

    """

    dg = dataset
    dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

    M_dict = {}

    # FOR EACH band:
    for band in ["J", "H", "K"]:
        mask = dat[f"{band}APERMAG3"].mask
        times = dat["MEANMJDOBS"][~mask]
        mags = dat[f"{band}APERMAG3"][~mask]
        errs = dat[f"{band}APERMAG3ERR"][~mask]

        try:
            M_dict[band] = compute_M_band(times, mags, errs)
        except IndexError:
            M_dict[band] = np.nan

    return M_dict


def compute_M_band(times, mags, errs):

    p_10 = np.percentile(mags, 10)
    p_90 = np.percentile(mags, 90)

    m_10 = np.mean(mags[(mags < p_10) | (mags > p_90)])
    m_med = np.median(mags)

    sigma_m = np.std(mags)

    M = (m_10 - m_med) / (sigma_m)

    return M
