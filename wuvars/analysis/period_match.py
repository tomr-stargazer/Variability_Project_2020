"""
A simple script which pulls in a handful of literature period surveys and 
compares our periods to them.

References the IPython notebook `Exploring literature periods.ipynb` a lot.

30 Sep 2024

"""

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
f1_ra_deg = f1[
    "RA\n            .",
]
f1_de_deg = f1[
    "Dec.\n            .",
]

f1_periods = f1["P\n            ."]

f1_litpertable = LitPeriodTable(
    f1,
    f1_ra_deg,
    f1_de_deg,
    f1_periods,
    short_name="F16",
    long_name="Fritzewski et al. 2016",
)


# All of this prep is somewhat sloppy.


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
        (f1_litpertable, w1_litpertable, g1_litpertable),
    ]

    symbols = ['o', 'o', 'o', 'v']
    ms_list = [3.5, 3, 2.5]

    for _periods, _coords, _name, litpertable_list in zip(
        cluster_periods, cluster_coords, cluster_names, cluster_litpertable_lists
    ):
        print(f"{_name}")
        print("")

        fig = plt.figure(figsize=(5, 5), dpi=150)
        ax = fig.add_subplot(111)
        ax.set_aspect("equal")

        period_plots = []

        for litpertable, symbol, ms in zip(litpertable_list, symbols, ms_list):

            idx, d2d, d3d = _coords.match_to_catalog_sky(litpertable.coords)

            max_sep = 1.0 * u.arcsec
            sep_constraint = d2d < max_sep

            print(f"Matches to {litpertable.long_name}: {np.sum(sep_constraint)}")
            print(f"Closest match: {np.min(d2d).to(u.arcsec):.2f}")
            print(
                f"Farthest used match: {np.max(d2d[sep_constraint]).to(u.arcsec):.2f}"
            )

            n_overlap = np.sum(
                np.isfinite(_periods["Period"][sep_constraint])
                & np.isfinite(litpertable.period[idx[sep_constraint]])
            )
            print(f"Number of overlapping periods: {n_overlap}")

            print("")

            matches = litpertable.table[idx[sep_constraint]]

            period_plots.extend(
                ax.plot(
                    _periods["Period"][sep_constraint],
                    litpertable.period[idx[sep_constraint]],
                    marker=symbol,
                    linestyle='None',
                    ms=ms,
                    label=litpertable.long_name + f" (n={n_overlap})",
                    zorder=20,
                )
            )

        # Create a legend for the first line.
        first_legend = ax.legend(handles=period_plots, loc="upper left")

        # Add the legend manually to the Axes.
        ax.add_artist(first_legend)

        alias_plots = []

        alias_plots.extend(
            ax.plot([0.01, 100], [0.01, 100], "k--", lw=0.25, label="$y=x$",)
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
        alias_plots.extend(ax.plot(xs_b, 1 / ys2_inv, "k:", lw=0.1,))
        alias_plots.extend(ax.plot(xs_c, 1 / ys3_inv, "k:", lw=0.1,))

        ax.axvline(1, color="k", linestyle=":", lw=0.25)
        ax.axhline(1, color="k", linestyle=":", lw=0.25)

        ax.set_xlim(0, 17)
        ax.set_ylim(0, 17)

        ax.set_xlabel("Our period (d)")
        ax.set_ylabel("Literature period (d)")

        ax.set_title(f"Literature period comparisons in {_name}")

        ax.legend(handles=alias_plots, loc="lower right")
