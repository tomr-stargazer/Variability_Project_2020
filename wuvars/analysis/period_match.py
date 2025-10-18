"""
A simple script which pulls in a handful of literature period surveys and
compares our periods to them.

References the IPython notebook `Exploring literature periods.ipynb` a lot.

30 Sep 2024

"""

import pdb

import astropy.table
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from wuvars.analysis.load_periodics_v4 import (ic_periods, ngc_periods,
                                               select_periodic_variables_v4)


class LitPeriodTable(object):
    def __init__(
        self, table, ra, dec, period, coords=None, short_name=None, long_name=None
    ):

        self.table = table
        self.ra = ra
        self.dec = dec

        if coords is None:
            coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg)).flatten()

        self.coords = coords
        self.period = period

        self.short_name = short_name
        self.long_name = long_name


# 1. Rebull 2015
r1 = astropy.table.Table.read(
    "../literature_periods/Rebull2015_aj520633t1_mrt.txt", format="ascii"
)

r1_periods = np.nanmean(
    np.ma.stack([r1["Period-i1"], r1["Period-i2"], r1["Period-i1i2"]]), axis=0
)
r1_litpertable = LitPeriodTable(
    r1,
    r1["RAdeg"],
    r1["DEdeg"],
    r1_periods,
    short_name="R15",
    long_name="Rebull et al. 2015",
)

# 2. Ghosh et al 2021
g1 = astropy.table.Table.read(
    "../literature_periods/Ghosh21_Fast photometric variability of very low mass stars in IC 348_ detection of superflare in an M dwarf _ Oxford Academic.html",
    format="ascii.html",
)

g1_ra_str = g1["RA\n            ."]
g1_de_str = g1["Dec.\n            ."]

g1_c = SkyCoord(ra=g1_ra_str, dec=g1_de_str, unit=(u.hourangle, u.deg))

g1_per_str = [
    (x[:-1]) if type(x) == np.str_ else "" for x in g1["Period\n            ."]
]
g1_per_str2 = ["nan" if x == "" else x for x in g1_per_str]
g1_per_hr = [float(x) if x != "–" else np.nan for x in g1["Period (h)\n            ."]]

g1_periods = np.nanmax(
    [np.array(g1_per_str2).astype(float), np.array(g1_per_hr) / 24], axis=0
)

g1_litpertable = LitPeriodTable(
    g1,
    g1_c.ra,
    g1_c.dec,
    g1_periods,
    coords=g1_c,
    short_name="G21",
    long_name="Ghosh et al. 2021",
)

# 3. Wang et al 2023
w1 = astropy.table.Table.read(
    "../literature_periods/Wang2023_raaacd58bt3_ascii.txt",
    format="ascii",
    delimiter="\t",
    header_start=2,
    data_start=4,
    data_end=238,
)

w1_ra_deg = w1["R.A."].astype(float)
w1_de_deg = w1["Decl."].astype(float)

w1_p = (w1["Category"] == "P") | (w1["Category"] == "P ^c")
w1_periods = w1["Timescale"]
w1_periods[~w1_p] = np.nan

w1_litpertable = LitPeriodTable(
    w1,
    w1_ra_deg,
    w1_de_deg,
    w1_periods,
    short_name="W23",
    long_name="Wang et al. 2023",
)

# 4. Fritzewski et al. 2016
f1 = astropy.table.Table.read(
    "../literature_periods/Fritzewski16_Long-term photometry of IC 348 with the Young Exoplanet Transit Initiative network _ Oxford Academic.html",
    format="ascii.html",
)
f1_ra_deg = f1["RA\n            .",]
f1_de_deg = f1["Dec.\n            .",]

f1_periods = f1["P\n            ."]

f1_litpertable = LitPeriodTable(
    f1,
    f1_ra_deg,
    f1_de_deg,
    f1_periods,
    short_name="F16",
    long_name="Fritzewski et al. 2016",
)

fla_1 = astropy.table.Table.read(
    "../literature_periods/Flaherty13_aj456236t1_mrt.txt", format="ascii"
)

fla_4 = astropy.table.Table.read(
    "../literature_periods/Flaherty13_aj456236t4_ascii.txt",
    format="ascii",
    delimiter="\t",
    header_start=2,
    data_start=4,
    data_end=25 + 6,
)

fla4_pers = []

for row in fla_4:

    sigs = [row["[3.6] Significance"], row["[4.5] Significance"]]
    choice = np.argmax(sigs)
    per = [row["[3.6] Period"], row["[4.5] Period"]][choice]

    fla4_pers.append(per)

fla_4_lrll = np.array([(x[5:]) for x in fla_4["Star ID"]])
fla1_pers = []

for row in fla_1:
    if "LRLL " + str(row["LRLL"]) in fla_4["Star ID"]:
        fla_4_match_index = np.int(np.where(row["LRLL"] == fla_4_lrll)[0])
        fla1_pers.append(fla4_pers[fla_4_match_index].rstrip("^a"))
    else:
        fla1_pers.append(np.nan)


fla1_litpertable = LitPeriodTable(
    fla_1,
    fla_1["RAdeg"],
    fla_1["DEdeg"],
    np.array(fla1_pers).astype(np.float),
    short_name="F13",
    long_name="Flaherty et al. 2013",
)

if __name__ == "__main__":

    # testing this out on one case

    ic_coords = SkyCoord(ra=ic_periods["RA_deg"], dec=ic_periods["DE_deg"], unit=u.deg)
    ngc_coords = SkyCoord(
        ra=ngc_periods["RA_deg"], dec=ngc_periods["DE_deg"], unit=u.deg
    )

    # start with NGC

    cluster_coords = [ngc_coords, ic_coords]
    cluster_periods = [ngc_periods, ic_periods]
    cluster_names = ["NGC 1333", "IC 348"]

    cluster_litpertable_lists = [
        (r1_litpertable, w1_litpertable),
        (f1_litpertable, w1_litpertable, fla1_litpertable, g1_litpertable),
    ]

    symbols = ["o", "o", "o", "o"]
    ms_list = [3.5, 3, 2.5, 2]

    for _periods, _coords, _name, litpertable_list in zip(
        cluster_periods, cluster_coords, cluster_names, cluster_litpertable_lists
    ):
        print(f"{_name}")
        print(f"{np.sum(np.isfinite(_periods['Period']))} periodic objects")
        print("")

        fig = plt.figure(figsize=(5, 5), dpi=150)
        ax = fig.add_subplot(111)
        ax.set_aspect("equal")

        period_plots = []

        for litpertable, symbol, ms in zip(litpertable_list, symbols, ms_list):

            idx, d2d, d3d = _coords.match_to_catalog_sky(litpertable.coords)

            max_sep = 1.0 * u.arcsec
            sep_constraint = d2d < max_sep

            litpertable.idx = idx
            litpertable.sep_constraint = sep_constraint

            print(f"  Matches to {litpertable.long_name}: {np.sum(sep_constraint)}")
            print(
                f"    Of these, {np.sum(np.isfinite(litpertable.period[idx[sep_constraint]]))} were identified by {litpertable.short_name} as periodic"
            )
            print(f"  Closest match: {np.min(d2d).to(u.arcsec):.2f}")
            print(
                f"  Farthest used match: {np.max(d2d[sep_constraint]).to(u.arcsec):.2f}"
            )

            litpertable.overlaps = (
                np.isfinite(_periods["Period"][sep_constraint])
                & np.isfinite(litpertable.period[idx[sep_constraint]])
            ) == True

            litpertable.alt_overlaps = (
                np.isfinite(_periods["Period"])
                & np.isfinite(litpertable.period[idx])
                & sep_constraint
            ) == True

            litpertable.was_periodic = (
                np.isfinite(litpertable.period[idx]) & sep_constraint
            ) == True

            n_overlap = np.sum(litpertable.overlaps)
            print(f"  Number of overlapping periods: {n_overlap}")

            print("")

            # Not sure I'm calculating this right.
            # What exactly are we measuring?

            # There are two things we might measure.

            # One is the literal "how many periodic objects, total, did I find
            # that are not literally listed as periodic in their catalog?"
            # and this would be a sum of the following two groups:
            #   (a) our-periodics that positionally matched to a non-periodic in their catalog;
            #   (b) our-periodics that did not match to their catalog.

            # It should literally be # of our found periods minus the number of periodic objects that matched to a periodic object

            # Break it down a little more finely...
            # a row in our table is flagged as a new period IF
            #   (a) it has a finite period in our data AND EITHER
            #      (b) it is not in the overlap period set among positional matches, OR
            #      (c) it did not match to an object in their catalog at all

            # I think I'm running into a challenge re: indexing on that latter bit.
            # Because I want everything indexed to
            our_new_periods = np.isfinite(_periods["Period"]) & (
                ~litpertable.alt_overlaps | ~sep_constraint
            )
            assert np.sum(our_new_periods) == np.sum(
                np.isfinite(_periods["Period"])
            ) - np.sum(litpertable.overlaps)

            our_new_periods_among_matches = np.isfinite(
                _periods["Period"][sep_constraint]
            ) & (~litpertable.overlaps)

            periods_we_missed = ~np.isfinite(
                _periods["Period"][sep_constraint]
            ) & np.isfinite(litpertable.period[idx[sep_constraint]])

            print(
                f"  Number of our new periods [relative to {litpertable.long_name}]: {np.sum(our_new_periods)}"
            )
            print(
                f"    Of which {np.sum(our_new_periods_among_matches)} were affirmatively not identified by {litpertable.short_name} as periodic"
            )
            print(
                f"  Number of periods we missed [relative to {litpertable.long_name}]: {np.sum(periods_we_missed)}"
            )
            print("")
            print("")

            matches = litpertable.table[idx[sep_constraint]]
            litpertable.matches = matches

            # pdb.set_trace()

            period_plots.extend(
                ax.plot(
                    _periods["Period"][sep_constraint],
                    litpertable.period[idx[sep_constraint]],
                    marker=symbol,
                    linestyle="None",
                    ms=ms,
                    label=litpertable.long_name + f" (n={n_overlap})",
                    zorder=20,
                )
            )

        ###
        # Here's where it makes sense to do some math on the periodic recovery fraction / number of new periodics.

        # total number of sources
        print(f"Total number of sources (periodic and not) in {_name}: {len(_periods)}")

        is_periodic = np.isfinite(_periods["Period"])
        n_is_periodic = np.sum(is_periodic)
        print(f"Number of periodics: {n_is_periodic}")

        overlaps_any = np.vstack(
            [litpertable.alt_overlaps for litpertable in litpertable_list]
        )
        overlaps = np.any(overlaps_any, axis=0)

        # was periodic in at least one of the literatures
        # may require some kind of list comprehension over listpertable_list

        was_periodic_any = np.vstack(
            [litpertable.was_periodic for litpertable in litpertable_list]
        )
        was_periodic = np.any(was_periodic_any, axis=0)

        print(
            f"  Of these, {np.sum(overlaps)} were previously identified as periodic by {[x.short_name for x in litpertable_list]}"
        )

        new_periodics = np.isfinite(_periods["Period"]) & ~was_periodic
        print("******** ********** ******** ****** ********")
        print(f"Number of **new** periodics we found: {np.sum(new_periodics)}")
        print("******** ********** ******** ****** ********")

        missed_periodics = ~np.isfinite(_periods["Period"]) & was_periodic

        print("")
        print(f"Number of known periodics we missed: {np.sum(missed_periodics)}")

        print("")
        print("")

        # pdb.set_trace()

        # Create a legend for the first line.
        first_legend = ax.legend(handles=period_plots, loc="upper left")

        # Add the legend manually to the Axes.
        ax.add_artist(first_legend)

        alias_plots = []

        alias_plots.extend(
            ax.plot(
                [0.01, 100],
                [0.01, 100],
                "k--",
                lw=0.25,
                label="$y=x$",
            )
        )
        alias_plots.extend(
            ax.plot(
                [0.01, 100],
                [0.01 * 2, 100 * 2],
                "k--",
                lw=0.1,
                label="$y=2x$, $y=x/2$",
            )
        )
        alias_plots.extend(ax.plot([0.01, 100], [0.01 / 2, 100 / 2], "k--", lw=0.1))

        # do the hyperbolic ones
        xs_a = np.linspace(1.01, 20, 200)
        xs_b = np.linspace(0.01, 20, 200)
        xs_c = np.linspace(0.01, 0.99, 50)
        ys1_inv = 1 - 1 / xs_a
        ys2_inv = 1 + 1 / xs_b
        ys3_inv = 1 / xs_c - 1

        alias_plots.extend(
            ax.plot(
                xs_a,
                1 / ys1_inv,
                "k:",
                lw=0.1,
                label=r"$y = \rm{abs} \left( \pm \left( \frac{1}{x} \pm 1 \right) \right ) ^{-1}$",
            )
        )
        alias_plots.extend(
            ax.plot(
                xs_b,
                1 / ys2_inv,
                "k:",
                lw=0.1,
            )
        )
        alias_plots.extend(
            ax.plot(
                xs_c,
                1 / ys3_inv,
                "k:",
                lw=0.1,
            )
        )

        ax.axvline(1, color="k", linestyle=":", lw=0.25)
        ax.axhline(1, color="k", linestyle=":", lw=0.25)

        ax.set_xlim(0, 17)
        ax.set_ylim(0, 17)

        ax.set_xlabel("Our period (d)")
        ax.set_ylabel("Literature period (d)")

        ax.set_title(f"Literature period comparisons in {_name}")

        ax.legend(handles=alias_plots, loc="lower right")
