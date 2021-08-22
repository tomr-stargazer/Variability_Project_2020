"""
Period recovery: a module dealing with the injection of artificial periodic signals
into nonvariable light curves and attempting to retrieve the injected signals using 
our period-finding tools.

A good way to characterize the variability that we see.

I'm imagining a grid like the following:

amplitudes:
- 0.005 mag (below ~half our noise limit)
- 0.01 mag
- 0.02 mag
- 0.05 mag
- 0.1 mag
- 0.5 mag

periods:
- 1 hour
- 6 hour
- 12 hour
- close to 1 day but not exactly
- 2 d
- 10 d
- 50 d
- 100 d

substrate star:
- ~5 random stars in each brightness bin, spaced such that all 4 "footprints" are covered
- brightness bins where the photometric noise is "optimal", "2x optimal", "5x optimal", "10x optimal"

then plot (input period, output period) for each of the above. Look at periodograms etc.

Do that for each SFR.

"""

import numpy as np
import pdb

# import matplotlib.pyplot as plt
# import pandas as pd
from astropy.timeseries import LombScargle

from wuvars.analysis.periods import freq


def recovery_score(dg, sids, periods, amplitudes):

    scores = []
    found_periods = []
    faps = []

    for period in periods:

        _scores = []
        _found_periods = []
        _faps = []

        for amp in amplitudes:

            __found_periods = []
            __faps = []
            correct = 0

            for sid in sids:

                dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

                mask = dat["KAPERMAG3"].mask
                times = dat["MEANMJDOBS"][~mask]
                mags = dat["KAPERMAG3"][~mask]
                errs = dat["KAPERMAG3ERR"][~mask]

                sin_mags = amp * np.sin(2 * np.pi / period * times) + mags

                ls = LombScargle(times, sin_mags, dy=errs)
                power = ls.power(freq)

                min_freq = 1 / 100
                power[freq < min_freq] = 0

                fmax = freq[np.nanargmax(power)]
                fap = ls.false_alarm_probability(np.nanmax(power))

                found_period = 1 / fmax

                if np.abs(found_period - period) / period < 0.01:
                    correct += 1
                    print(
                        f"A={amp:.2f} mag. Correct period: {found_period:.2f} v. {period:.2f}"
                    )
                else:
                    print(
                        f"A={amp:.2f} mag. Incorrect period: {found_period:.2f} v. {period:.2f}"
                    )

                # pdb.set_trace()

                __found_periods.append(found_period)
                __faps.append(fap)

            score = correct / len(sids)
            _scores.append(score)
            _found_periods.append(__found_periods)
            _faps.append(__faps)

        scores.append(_scores)
        found_periods.append(_found_periods)
        faps.append(_faps)

    return scores, found_periods, faps
