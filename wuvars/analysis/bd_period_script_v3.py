"""
It's just a script. Loads up the relevant data, does periodograms, 
saves those periodograms as images, and saves the periodogram best-fit output
to the table columns.

"""

import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import pdb

import astropy.table
from astropy.timeseries import LombScargle

from wuvars.data import spreadsheet, photometry, quality_classes
from wuvars.analysis.periods import f_max, f_min, N_eval
from wuvars.analysis.alias_hunting import find_aliases, find_m_aliases, find_n_aliases
from wuvars.plotting.lightcurve import simple_phased_lc_scatter_gridspec
from wuvars.plotting.lightcurve_helpers import orion_cmap

plt.style.use("seaborn-whitegrid")
warnings.filterwarnings("ignore")

freq = np.linspace(f_min, f_max, N_eval)

onc = {}
ngc = {}
ic = {}
dicts = [onc, ngc, ic]

onc["dat"] = photometry.group_wserv_v2(photometry.load_wserv_v2(5))
onc["q"] = quality_classes.load_q(5)
onc["spread"] = spreadsheet.load_wserv_v2(5)
onc["cmap"] = orion_cmap

ngc["dat"] = photometry.group_wserv_v2(photometry.load_wserv_v2(7))
ngc["q"] = quality_classes.load_q(7)
ngc["spread"] = spreadsheet.load_wserv_v2(7)
ngc["cmap"] = "jet"

ic["dat"] = photometry.group_wserv_v2(photometry.load_wserv_v2(8))
ic["q"] = quality_classes.load_q(8)
ic["spread"] = spreadsheet.load_wserv_v2(8)
ic["cmap"] = "jet_r"

region_names = ["ONC", "NGC1333", "IC348"]
short_names = ["ONC", "NGC", "IC"]

bands = ["J", "H", "K"]


def write_periods_and_generate_periodograms():
    for dict_, name, short_name in zip(dicts, region_names, short_names):

        dg = dict_["dat"]

        # don't load the corresponding table. We'll combine them later.
        dir_ = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes"
        source_properties = astropy.table.Table.read(
            os.path.join(dir_, short_name + "_source_properties.csv")
        )

        # place to put periodograms
        results_dir = (
            "/Users/tsrice/Documents/Variability_Project_2020/Results/periodograms"
        )

        for i, sid in enumerate(source_properties["SOURCEID"]):
            dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

            fig, axes = plt.subplots(
                nrows=3, ncols=1, sharex=True, figsize=(10, 10), dpi=200
            )

            for (band, ax) in zip(bands, axes):

                mask = dat[f"{band}APERMAG3"].mask
                times = dat["MEANMJDOBS"][~mask]
                mags = dat[f"{band}APERMAG3"][~mask]
                errs = dat[f"{band}APERMAG3ERR"][~mask]

                # compute a lomb scargle
                try:
                    ls = LombScargle(times, mags, dy=errs)
                    power = ls.power(freq, assume_regular_frequency=True).value

                    min_freq = 1 / 100
                    power[freq < min_freq] = 0
                    power[np.abs(freq-1) < 0.01] = 0

                    powermax = np.nanmax(power)
                    fmax = freq[np.nanargmax(power)]
                    fap = ls.false_alarm_probability(np.nanmax(power)).value

                    amp = np.sqrt(np.sum(ls.model_parameters(fmax)[1:3]**2))
                    # print(ls.model_parameters(fmax), amp)

                    # find aliases
                    m_aliases = np.array(find_m_aliases(fmax))
                    n_aliases = np.array(find_n_aliases(fmax))

                    # plot it in period space.
                    ax.plot(1 / freq, power, "k", lw=0.8, rasterized=True)
                    ax.set_ylabel(f"{band} power")
                    if band == "K":
                        ax.set_xlabel("Period (d)")
                        ax.set_xscale("log")
                        ax.set_xlim(1 / 24, 100)
                    # ax.axvline(1, ls="--", lw=0.5, alpha=0.3)

                    # tag aliases
                    ax.set_ylim(0, powermax)
                    ymin, ymax = ax.get_ylim()

                    outside_line = ax.vlines(
                        1 / fmax, ymax * 1.05, ymax * 1.1, color="k"
                    )
                    red_line = ax.vlines(
                        1 / m_aliases, ymax * 1.0, ymax * 1.05, color="C3"
                    )
                    blue_line = ax.vlines(
                        1 / n_aliases, ymax * 1.01, ymax * 1.06, color="C9"
                    )

                    outside_line.set_clip_on(False)
                    red_line.set_clip_on(False)
                    blue_line.set_clip_on(False)

                except ValueError:
                    fmax = 1
                    amp = 0
                    fap = 1

                # update the table
                source_properties[f"period_{band}"][i] = 1 / fmax
                source_properties[f"per_amp_{band}"][i] = amp
                source_properties[f"per_fap_{band}"][i] = fap

            fig.suptitle(f"{name} {i:03d} SID: {sid}")

            min_fap = min([source_properties[f"per_fap_{band}"][i] for band in bands])
            # pdb.set_trace()

            if -np.log10(min_fap) > 5:
                # pdb.set_trace()
                print(f"{name} {i:03d} is a periodic candidate (FAP: {min_fap:.2e})")

                fig.savefig(
                    os.path.join(results_dir, name, f"{i:03d}_periodogram.png"),
                    bbox_inches="tight",
                )
            else:
                print(f"{name} {i:03d} is not a periodic candidate (FAP={min_fap:.2e})")
            plt.close(fig)

        print("Writing updated source_properties to:")
        print(os.path.join(dir_, short_name + "_source_properties_periods.csv"))
        source_properties.write(
            os.path.join(dir_, short_name + "_source_properties_periods.csv")
        )
    return None


def generate_folded_lightcurves():
    for dict_, name, short_name in zip(dicts, region_names, short_names):
        plt.style.use("default")
        dg = dict_["dat"]

        # don't load the corresponding table. We'll combine them later.
        dir_ = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes"
        source_properties = astropy.table.Table.read(
            os.path.join(dir_, short_name + "_source_properties_periods.csv")
        )

        # place to put periodograms
        results_dir = "/Users/tsrice/Documents/Variability_Project_2020/Results/folded_lightcurves"

        for i, sid in enumerate(source_properties["SOURCEID"]):

            min_fap = min([source_properties[f"per_fap_{band}"][i] for band in bands])
            if -np.log10(min_fap) > 5:
                print(f"{name} {i:03d} is a periodic candidate (FAP={min_fap:.2e})")

                dat = dg.groups[dg.groups.keys["SOURCEID"] == sid]

                for band in bands:

                    period = source_properties[f"period_{band}"][i]
                    fap = source_properties[f"per_fap_{band}"][i]

                    # make its 3 folded lightcurves, colored by time
                    fig = simple_phased_lc_scatter_gridspec(
                        dat, sid, period, cmap=dict_["cmap"]
                    )
                    fig.suptitle(
                        f"{band} period: {period:.2f}d | {name} {i:03d} SID: {sid}"
                    )

                    fig.savefig(
                        os.path.join(
                            results_dir, name, f"{i:03d}_{band}_folded_lc.png"
                        ),
                        bbox_inches="tight",
                    )
                    plt.close(fig)
            else:
                print(f"{name} {i:03d} is not a periodic candidate (FAP={min_fap:.2e})")

    return None


def append_images():
    from subprocess import Popen

    for dict_, name, short_name in zip(dicts, region_names, short_names):

        dir_ = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes"
        source_properties = astropy.table.Table.read(
            os.path.join(dir_, short_name + "_source_properties_periods.csv")
        )

        # place to put periodograms
        lc_dir = os.path.join(
            "/Users/tsrice/Documents/Variability_Project_2020/Results/folded_lightcurves",
            name,
        )
        pgram_dir = os.path.join(
            "/Users/tsrice/Documents/Variability_Project_2020/Results/periodograms",
            name,
        )

        for i, sid in enumerate(source_properties["SOURCEID"]):

            min_fap = min([source_properties[f"per_fap_{band}"][i] for band in bands])
            if -np.log10(min_fap) > 5:

                cmd_str1 = f"convert {lc_dir}/{i:03d}_J_folded_lc.png {lc_dir}/{i:03d}_H_folded_lc.png {lc_dir}/{i:03d}_K_folded_lc.png -append {lc_dir}/{i:03d}_JHK_folded_lc.png"
                cmd_str2 = f"convert {lc_dir}/{i:03d}_JHK_folded_lc.png \\( {pgram_dir}/{i:03d}_periodogram.png -resize 78% \\) -gravity center +append {pgram_dir}/{i:03d}_JHK_lc_pgram_composite.png"

                # print(cmd_str1)
                # print(cmd_str2)

                Popen(cmd_str1+'\n'+cmd_str2, shell=True)
                # Popen(cmd_str2, shell=True)


if __name__ == "__main__":
    if True:
        write_periods_and_generate_periodograms()
        generate_folded_lightcurves()
        append_images()
