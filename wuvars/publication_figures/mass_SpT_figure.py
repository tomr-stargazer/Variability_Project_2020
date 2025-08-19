"""
This is a script for making a mass-SpT figure for publication.

I had seriously prototyped such a figure here: 

/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/Prototyping a mass-SpT figure.ipynb

and also, in part, here:

/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/Playing%20with%20Baraffe%20Isochrones.ipynb

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from wuvars.analysis.spectral_type_to_number import (get_num_from_SpT,
                                                     get_SpT_from_num)
from wuvars.analysis.spectral_type_to_temperature import (
    get_SpT_from_Teff, get_SpT_from_Teff_HH14)
from wuvars.publication_figures.HR_diagram_NGC_IC import (
    extract_age_array_Gyr, load_isochrone_generic)

mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.right"] = True

csv_filename = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/publication_figures/TeffSpT_relations.csv"
df = pd.read_csv(csv_filename)


def make_mass_SpT_figure():

    ages = [1, 2, 3, 5]

    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    plt.grid(True, linestyle="--", lw=0.25)

    lw_array = [0.5, 0.75, 1.0, 1.25]
    ms_array = [2, 3, 4, 5]
    colors = [f"C{x}" for x in range(6)]

    for age, lw, ms in list(zip(ages, lw_array, ms_array))[::-1]:
        iso = load_isochrone_generic(age)
        mass = iso[:, 0]
        teff = iso[:, 1]

        spt_array = df["SpTnum"]
        for col, color in zip(df.columns[2:-1], colors):
            #        print(col)

            nonnan = ~np.isnan(df[col].values)

            SpT_interp = lambda teff: np.interp(
                teff,
                df[col].values[nonnan][::-1],
                spt_array.values[nonnan][::-1],
                left=np.nan,
                right=np.nan,
            )

            SpTs = SpT_interp(teff)
            #         print(col, SpTs)

            small = mass <= 0.75
            ax.plot(
                SpTs[small], mass[small], ":", color=color, lw=lw, alpha=0.2, zorder=ms
            )
            ax.plot(
                SpTs[small],
                mass[small],
                ".",
                color=color,
                lw=lw,
                ms=ms,
                label=col,
                zorder=ms,
            )

        if age == ages[-1]:
            ax.legend()

    ax.axhline(0.08, zorder=-10, alpha=0.1, lw=5)

    ax.set_xlim(-0.25, 12.25)
    ax.minorticks_on()

    # gets all the xticks, converts them into proper spectral type strings, then puts them as labels
    ax.set_xticklabels(np.array([get_SpT_from_num(int(x)) for x in ax.get_xticks()]))

    ax.set_xlabel("Spectral Type")
    ax.set_ylabel("Mass (M$_{\odot}$)")
    ax.semilogy()

    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_ticks([0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8])

    # Now to workshop the baby axes

    ax_bounds = ax.get_position().bounds
    left = ax_bounds[0]

    left_buffer = left
    right_buffer = 1 - (ax_bounds[0] + ax_bounds[2])
    total = 1 - (left_buffer + right_buffer)
    sum_width = total - (left_buffer + right_buffer) / 2
    half_width = sum_width / 2

    bottom = -half_width
    width = half_width
    height = half_width

    ax_left = fig.add_axes([left, bottom, width, height])
    ax_right = fig.add_axes(
        [left + width + (left_buffer + right_buffer) / 2, bottom, width, height]
    )

    ax_left.tick_params(axis="both", labelsize=6)
    ax_right.tick_params(axis="both", labelsize=6)

    # LEFT

    ax_left.grid(True, linestyle="--", lw=0.25)

    for age, lw in zip(ages, lw_array):

        myr_x = load_isochrone_generic(age)

        mass = myr_x[:, 0]
        teff = myr_x[:, 1]
        Mk = myr_x[:, -1]
        M = myr_x[:, 0]

        u = M <= 0.8

        ax_left.plot(mass[u], teff[u], "k", lw=lw, label=f"{age} Myr")

    ax_left.set_ylabel(r"$T_{\rm{eff}}$ (K)", fontsize=8)
    ax_left.set_xlabel("Mass (M$_{\odot})$", fontsize=8)
    ax_left.legend(fontsize=6)

    ax_left.semilogx()
    ax_left.xaxis.set_major_formatter(ScalarFormatter())
    ax_left.xaxis.set_ticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5])
    ax_left.axvline(0.08, zorder=-10, alpha=0.1, lw=5)
    ax_left.set_xlim(0.01, 0.81)

    # RIGHT

    ax_right.grid(True, linestyle="--", lw=0.25)

    spt_array = df["SpTnum"]
    for col in df.columns[2:-1]:
        ax_right.plot(spt_array, df[col].values, ".-", lw=0.75, ms=5, label=col)

    ax_right.legend(fontsize=6)

    ax_right.set_xticks(ax.get_xticks())
    ax_right.set_xlim(ax.get_xlim())
    ax_right.set_xticklabels(
        np.array([get_SpT_from_num(int(x)) for x in ax_right.get_xticks()])
    )
    ax_right.set_ylim(1500, 4000)
    ax_right.set_ylabel(r"$T_{\rm{eff}}$ (K)", fontsize=8)
    ax_right.set_xlabel("Spectral Type", fontsize=8)

    plt.savefig("XXX_spectral_type_v_mass_logy_altered.pdf")

    return fig


if __name__ == "__main__":

    fig = make_mass_SpT_figure()
