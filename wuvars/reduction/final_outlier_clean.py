"""
This is a script that will attempt to remove the "final" outliers in 
graded-clipped-scrubbed-whatever data.

Here's what I wrote in `immediate post cool stars notes`:

> I’m going to resolve the issue in NGC 1333 of conspicuous “outlier” 
photometric datapoints. 
> I have a hunch that some combination of “grade” filtering plus 
“check if this data point is more than three times the photometric uncertainty 
from its neighbors plus if (e.g.) H and K are also not deviant” will handle it. 
I believe this is scientifically justifiable even if it is slightly ad hoc - 
it is, essentially, me formalizing the rule I already use to “by eye” judge whether 
there is any reliable physical meaning in any given “large outlier” data point. 
It is not simply “removing data that I don’t like”.
> I’ll also check to see if IC 348 needs a similar treatment.


"""

def calculate_diff(dat, sid, band="J", verbose=False):

    dt = dat.groups[dat.groups.keys["SOURCEID"] == sid]

    # dt

    j_vals = dt['JAPERMAG3']
    j_sigs = dt['JAPERMAG3ERR']

    diff = np.zeros_like(j_vals)

    l = len(j_vals)

    for i in range(l):

        if verbose:
            print("i=", i)

        left_diff = np.nan
        right_diff = np.nan
        diff[i] = np.nan

        if not np.ma.is_masked(j_vals[i]):

            a_l = 1

            while (i-a_l) >= 0:

                left_diff = (j_vals[i-a_l] - j_vals[i]) / j_sigs[i]

                if verbose:
                    print("Left diff:", left_diff, "a_l = ", a_l)

                if i==8:
                    pdb.set_trace()

                if not np.ma.is_masked(left_diff):
                    if verbose:
                        print("  Breaking -  we accept left_diff")
                    break
                else:
                    if verbose:
                        print("  Not breaking - we don't accept left_diff")
                    a_l += 1
                    if verbose:
                        print("  a_l = ", a_l)
            if verbose:
                print("a_l: ", a_l, "i-a_l=", i-a_l)
                print("")

            a_r = 1

            while (i+a_r) < l:

                right_diff = (j_vals[i+a_r] - j_vals[i]) / j_sigs[i]

                if not np.ma.is_masked(right_diff):
                    print("Breaking")
                    break
                else:
                    print("Not breaking")
                    a_r += 1     
            if verbose:
                print("a_r: ", a_r, "i+a_r", i+a_r, "l=", l)

                print(left_diff, right_diff)

            try:
                diffs = np.array([left_diff, right_diff])
                abs_diffs = np.array([np.abs(left_diff), np.abs(right_diff)])
                diff_abs = min(abs_diffs)
                diff_sgn = np.sign( diffs[diff_abs==abs_diffs] )[0]

                diff[i] = diff_sgn * diff_abs
            except IndexError:
                diff[i] = np.nan
        else:
            if verbose:
                print(f"i={i}: J={j_vals[i]} is masked")

    #         if i>3:
    #             break
    return diff


def compute_nearness_metric(df):

    

def clean(df):
    """
    pseudocode

    for each star
        for each timestamp:

            for each band:

                ask:
                    * is the data point more than 3sigma away from both its neighbors on each side?
                    * is its grade below 95% (or some other threshold)?
                    * do the other two bands haev data "within" 3 sigma of their neighbors?


    so there's probably a way in which I could code up some metric all in one pass for 
    the dataset which essentially assesses the first one -
    in other words, "compute how far in sigmas this data point is from its nearest 
    neighbors on either side."

    """
