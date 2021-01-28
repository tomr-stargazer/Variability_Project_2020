"""
In the paper "The UKIRT wide field camera ZYJHK photometric system: calibration from 
2MASS" by Hodgkin et al. (2009), the authors derive an ad-hoc "correction" to the 
pipeline-generated photometric uncertainty error bars (which they note* are explicitly 
not meant to take calibration uncertainties into account), which looks like this:

  "M^2 = c E^2 + s^2

  where the systematic component s = 0.021 ± 0.001, 
  the constant of proportionality c = 1.082 ± 0.014. 

  In principle, these can be used to update the default pipeline error estimates..."

In our surveys, we have enough data to re-derive this kind of analysis (because we have
SO many repeated observations), so I'm going to do that. Especially because I've 
rejected the worst data using my cleaning program, I think there's an explicit value in
doing this step and KNOWING what the errors are. I could even say something like:

  Hodgkin et al. (2009) demonstrated the feasibility of an approach which adjusts the 
  pipeline-estimated error bars by empirically measured factors to produce more 
  realistic error bars. Here, we reproduce their analysis using our own data (after 
  the aforementioned cleaning procedures) to gain a more

* "The pipeline errors are calculated from the source counts and local background only, 
and do not include corrections for bad pixels and other detector artefacts. In 
addition, they deliberately do not include any contribution from the calibration 
procedure itself (i.e. the inherent uncertainties in the calibrators) or additional 
corrections for residual systematic effects in the calibration." - Section 3.1, p. 678

"""

from collections import defaultdict
import numpy as np

# we're gonna want to go through - at first glance - each magnitude "bin" and get a
# representative 'observed rms' value, as well as a representative 'estimated error' value.

# imagine this is WSERV11, and we've selected (as best we can) a pristine, non-variable
# sample. For each band, in each magnitude range, what is a representative
# observed rms? (and estimated rms?)

# ok, let's do this.


def decorrect_error(M, s=0.021, c=1.082):
    """
    M^2 = cE^2 + s^2 
    therefore
    E = sqrt( (M^2 - s^2) / c )
    """

    E = np.sqrt((M ** 2 - s ** 2) / c)

    return E


def bin_the_median_errorbars_by_magnitude(ds, bins=10, range=None, decorrect=True):
    """
    Takes a properly sanitized and de-variabled summary spreadsheet `ds`, 
    returns three arrays per band:
    - bin centers (magnitude array)
    - typical MEASURED rms error (rms array)
    - typical ESTIMATED rms error

    kwarg `decorrect`: use True if your input data originally had a Hodgkin 2009 
    applied to it (likely if you are me, Tom!). Use False if you went back and undid 
    that.

    """

    bands = ["J", "H", "K"]
    output_dict = defaultdict(dict)

    for b in bands:

        # create bin edges for the magnitude array
        mags = ds["median"][f"{b}APERMAG3"]

        if range is None:
            hist_range = np.nanmin(mags), np.nanmax(mags)
        else:
            hist_range = range

        hist, bin_edges = np.histogram(mags, bins=bins, range=hist_range)

        bin_lefts, bin_rights = bin_edges[:-1], bin_edges[1:]
        bin_centers = (bin_lefts + bin_rights) / 2

        rms = ds["std"][f"{b}APERMAG3"]
        err = ds["median"][f"{b}APERMAG3ERR"]
        if decorrect:
            err = decorrect_error(err)

        rms_list = []
        err_list = []

        for left, right, center in zip(bin_lefts, bin_rights, bin_centers):

            median_observed_rms_in_this_mag_bin = np.nanmedian(
                rms[(mags > left) & (mags < right)]
            )

            median_estimated_error_in_this_mag_bin = np.nanmedian(
                err[(mags > left) & (mags < right)]
            )

            rms_list.append(median_observed_rms_in_this_mag_bin)
            err_list.append(median_estimated_error_in_this_mag_bin)

        output_dict[b]["mag_bin_center"] = bin_centers
        output_dict[b]["measured_error"] = rms_list
        output_dict[b]["estimated_error"] = err_list

    return output_dict
