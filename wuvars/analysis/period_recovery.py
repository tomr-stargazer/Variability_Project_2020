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
from wuvars.analysis.alias_hunting import find_aliases


def recovery_score(dg, sids, periods, amplitudes):

    scores = []
    found_periods = []
    faps = []
    alias_scores = []

    for period in periods:

        _scores = []
        _found_periods = []
        _faps = []
        _alias_scores = []

        for amp in amplitudes:

            __found_periods = []
            __faps = []
            correct = 0
            alias_correct = 0

            for sid in sids:

                dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

                mask = dat["KAPERMAG3"].mask
                times = dat["MEANMJDOBS"][~mask]
                mags = dat["KAPERMAG3"][~mask]
                errs = dat["KAPERMAG3ERR"][~mask]

                sin_mags = amp * np.sin(2 * np.pi / period * times) + mags

                ls = LombScargle(times, sin_mags, dy=errs)
                power = ls.power(freq, assume_regular_frequency=True).value

                min_freq = 1 / 100
                power[freq < min_freq] = 0

                fmax = freq[np.nanargmax(power)]
                fap = ls.false_alarm_probability(np.nanmax(power)).value

                # pdb.set_trace()

                found_period = 1 / fmax

                if np.abs(found_period - period) / period < 0.02:
                    correct += 1
                    print(
                        f"A={amp:.2f} mag. Correct period: {found_period:.2f} v. {period:.2f}"
                    )
                else:
                    # CHECK ALIASES
                    aliases = find_aliases(1/period)

                    # LOOP THROUGH POSSIBLE ALIASES
                    for alias_freq in aliases:
                        alias_period = 1/alias_freq

                        if np.abs(found_period - alias_period) / alias_period < 0.02:
                            alias_correct += 1
                            print(
                                f"A={amp:.2f} mag. Aliased period: {found_period:.4f} v. {period:.4f}"
                            )
                            break
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
            alias_score = alias_correct / len(sids)
            _alias_scores.append(alias_score)

        scores.append(_scores)
        found_periods.append(_found_periods)
        faps.append(_faps)
        alias_scores.append(_alias_scores)

    return scores, found_periods, faps, alias_scores


def reassess_alias_score(spread, sids, periods, amplitudes, previous_output):
    """
    Assess how many of the periods land on aliases.

    """

    alias_scores = []
    found_periods = previous_output[1]
    faps = previous_output[2]

    for i, period in enumerate(periods):

        _new_scores = []
        _found_periods = found_periods[i]
        _faps = faps[i]

        for j, amp in enumerate(amplitudes):

            __found_periods = _found_periods[j]
            __faps = _faps[j]
            correct = 0
            skip = 0

            for k, sid in enumerate(sids):

                fap = __faps[k]
                found_period = __found_periods[k]

                if 0 < spread['variability']['Stetson_JHK'][sid] < 0.8:

                    # LOOP THROUGH POSSIBLE ALIASES
                    aliases = find_aliases(1/period)

                    for alias in aliases:
                        alias_period = 1/alias

                        if np.abs(found_period - alias_period) / alias_period < 0.01:
                            correct += 1
                        else:
                            pass
                    else:
                        skip +=1

                    # pdb.set_trace()

            new_score = correct / (len(sids) - skip)
            _new_scores.append(new_score)
        alias_scores.append(_new_scores)
    return alias_scores, found_periods, faps
