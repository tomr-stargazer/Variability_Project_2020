"""
In this script I am coding up my "HR diagram" figure.

We are going to do six panels:

two columns, one each for IC 348 and NGC 1333.

three rows:
- the HR diagram, aka K magnitude versus 
    I will show all objects included successfully in this study (so, 
    excluding non-positional matches and useless light curves).
    I will include age-appropriate 
- variability amplitude versus 


"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from wuvars.analysis.bd_matching_v3 import match_ic, match_ngc
from wuvars.analysis.spectral_type_to_number import (get_num_from_SpT,
                                                     get_SpT_from_num)
from wuvars.analysis.spectral_type_to_temperature import get_SpT_from_Teff_HH14
from wuvars.data import spreadsheet
from wuvars.publication_figures.make_figures import figure_export_path

# first - a proof of concept on the isochrones
filepath = (
    "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/BHAC15_iso.ukidss"
)
distmod = 5 * np.log10(300) - 5


def extract_age_array_Gyr(file):
    with open(file, "r") as fp:
        lines = fp.readlines()
        output = []
        for i, line in enumerate(lines):
            if "t (Gyr)" in line:
                line_age = line[14:-1]
                output.append(float(line_age))
    return np.array(output)


def find_line_given_age(age_Gyr, file):
    with open(file, "r") as fp:
        lines = fp.readlines()
        for i, line in enumerate(lines):
            if "t (Gyr)" in line:
                line_age = line[15:-1]
                if float(line_age) == age_Gyr:
                    #                     print("MATCH: ", line_age, age, i)
                    return i
    return None


def load_isochrone_generic(
    age_Myr,
    file="/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/BHAC15_iso.ukidss",
):

    start_line = find_line_given_age(age_Myr / 1e3, file)

    subtracted_rows = 0
    while subtracted_rows <= 30:
        try:
            myr_x = np.loadtxt(
                filepath, skiprows=start_line + 4, max_rows=30 - subtracted_rows
            )
            for i in range(subtracted_rows):
                pad = np.ones(11) * np.nan
                myr_x = np.vstack([pad, myr_x])
            return myr_x
        except ValueError:
            subtracted_rows += 1
        except TypeError:
            return None

    return myr_x


def mock_lalchand_figure():
    # Let's reproduce something very much like Lalchand et al. 2022
    # But - projected into the observable space for my data.

    fig, ax = plt.subplots(1, dpi=150)

    ages = [0.5, 2.0, 5.0, 10.0]

    for age in ages:

        myr_x = load_isochrone_generic(age)

        teff = myr_x[:, 1]
        Mk = myr_x[:, -1]
        M = myr_x[:, 0]

        u = M <= 0.3

        SpT = get_SpT_from_Teff_HH14(teff)

        plt.plot(SpT[u], Mk[u] + distmod, "k", lw=0.75, label=f"{age} Myr")

    all_ages_Gyr = extract_age_array_Gyr(filepath)

    masses = [0.3, 0.2, 0.1, 0.05]
    myr_0_5 = load_isochrone_generic(0.5)
    mass_array = myr_0_5[:, 0]

    for i, mass in enumerate(masses):
        # get the evolutionary track for all ages
        all_isochrones = []
        for age in all_ages_Gyr[all_ages_Gyr < 0.5]:
            all_isochrones.append(load_isochrone_generic(age * 1e3))
        myr_stack = np.vstack(all_isochrones).reshape(len(all_isochrones), 30, 11)

        mass_index = np.where(mass == mass_array)[0]
        mi = mass_index
        print(f"mass_index={mass_index}")

        teff_i = myr_stack[:, mi, 1]
        Mk_i = myr_stack[:, mi, -1]
        M_i = myr_stack[:, mi, 0]

        SpT_i = get_SpT_from_Teff_HH14(teff_i)

        plt.plot(SpT_i, Mk_i + distmod, "--", label=f"{mass} M$_{{\odot}}$")

    #     if M_i[0] <= 0.3:
    #         if M_i[0] in masslist:
    #             print("M = ", M_i[0])
    #             plt.plot(SpT_i, Mk_i + distmod, label=f"{M_i[0]} M$_{{\odot}}$")

    ax.set_xlim(2, 14)
    ax.invert_yaxis()
    ax.set_xlabel("Spectral Type (0=M0, 10=L0)")
    ax.set_ylabel("$K$ magnitude ($d=300$ pc)")
    plt.legend()

    return fig


lw_array = [0.25, 0.75, 1.5, 2]


def mock_lalchand_figure_v2():
    # Let's reproduce something very much like Lalchand et al. 2022
    # But - projected into the observable space for my data.

    fig, ax = plt.subplots(1, dpi=150)

    ages = [0.5, 2.0, 5.0, 10.0]

    for age, lw in zip(ages, lw_array):

        myr_x = load_isochrone_generic(age)

        teff = myr_x[:, 1]
        Mk = myr_x[:, -1]
        M = myr_x[:, 0]

        u = M <= 0.3

        SpT = get_SpT_from_Teff_HH14(teff)

        plt.plot(SpT[u], Mk[u] + distmod, "k", lw=lw, label=f"{age} Myr")

    all_ages_Gyr = extract_age_array_Gyr(filepath)

    masses = [0.2, 0.1, 0.05, 0.03, 0.01]
    myr_0_5 = load_isochrone_generic(0.5)
    mass_array = myr_0_5[:, 0]

    for i, mass in enumerate(masses):
        # get the evolutionary track for all ages
        all_isochrones = []
        for age in all_ages_Gyr[all_ages_Gyr < 0.5]:
            all_isochrones.append(load_isochrone_generic(age * 1e3))
        myr_stack = np.vstack(all_isochrones).reshape(len(all_isochrones), 30, 11)

        mass_index = np.where(mass == mass_array)[0]
        mi = mass_index
        print(f"mass_index={mass_index}")

        teff_i = myr_stack[:, mi, 1]
        Mk_i = myr_stack[:, mi, -1]
        M_i = myr_stack[:, mi, 0]

        SpT_i = get_SpT_from_Teff_HH14(teff_i)

        plt.plot(SpT_i, Mk_i + distmod, "--", label=f"{mass} M$_{{\odot}}$")

    #     if M_i[0] <= 0.3:
    #         if M_i[0] in masslist:
    #             print("M = ", M_i[0])
    #             plt.plot(SpT_i, Mk_i + distmod, label=f"{M_i[0]} M$_{{\odot}}$")

    ax.set_xlim(3, 13)
    ax.invert_yaxis()
    ax.set_xlabel("Spectral Type")

    # gets all the xticks, converts them into proper spectral type strings, then puts them as labels
    ax.set_xticklabels(np.array([get_SpT_from_num(int(x)) for x in ax.get_xticks()]))

    ax.set_ylabel("$K$ magnitude ($d=300$ pc)")
    plt.legend()

    return fig


def mock_lalchand_figure_v3():
    # Let's reproduce something very much like Lalchand et al. 2022
    # In Teff space!

    fig, ax = plt.subplots(1, dpi=150)

    ages = [0.5, 2.0, 5.0, 10.0]

    for age, lw in zip(ages, lw_array):

        myr_x = load_isochrone_generic(age)

        teff = myr_x[:, 1]
        Mk = myr_x[:, -1]
        M = myr_x[:, 0]

        u = M <= 0.3

        plt.plot(teff[u], Mk[u] + distmod, "k", lw=lw, label=f"{age} Myr")

    all_ages_Gyr = extract_age_array_Gyr(filepath)

    masses = [0.2, 0.1, 0.05, 0.03, 0.01]
    myr_0_5 = load_isochrone_generic(0.5)
    mass_array = myr_0_5[:, 0]

    for i, mass in enumerate(masses):
        # get the evolutionary track for all ages
        all_isochrones = []
        for age in all_ages_Gyr[all_ages_Gyr < 0.5]:
            all_isochrones.append(load_isochrone_generic(age * 1e3))
        myr_stack = np.vstack(all_isochrones).reshape(len(all_isochrones), 30, 11)

        mass_index = np.where(mass == mass_array)[0]
        mi = mass_index
        print(f"mass_index={mass_index}")

        teff_i = myr_stack[:, mi, 1]
        Mk_i = myr_stack[:, mi, -1]
        M_i = myr_stack[:, mi, 0]

        plt.plot(teff_i, Mk_i + distmod, "--", label=f"{mass} M$_{{\odot}}$")

    #     if M_i[0] <= 0.3:
    #         if M_i[0] in masslist:
    #             print("M = ", M_i[0])
    #             plt.plot(SpT_i, Mk_i + distmod, label=f"{M_i[0]} M$_{{\odot}}$")

    # ax.set_xlim(3, 13)
    ax.invert_yaxis()
    ax.invert_xaxis()
    ax.set_xlabel(r"$T_{\rm{eff}}$ (K)")
    ax.set_ylabel("$K$ magnitude ($d=300$ pc)")
    plt.legend()

    return fig


def mock_lalchand_figure_v4():
    # Let's reproduce something very much like Lalchand et al. 2022
    # Mass vs Teff!

    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)
    plt.grid(True, linestyle="--", lw=0.25)

    ages = [0.5, 2.0, 5.0, 10.0]

    for age, lw in zip(ages, lw_array):

        myr_x = load_isochrone_generic(age)

        mass = myr_x[:, 0]
        teff = myr_x[:, 1]
        Mk = myr_x[:, -1]
        M = myr_x[:, 0]

        u = M <= 0.8

        plt.plot(teff[u], mass[u], "k", lw=lw, label=f"{age} Myr")

    all_ages_Gyr = extract_age_array_Gyr(filepath)

    # masses = [0.2, 0.1, 0.05, 0.03, 0.01]
    # myr_0_5 = load_isochrone_generic(0.5)
    # mass_array = myr_0_5[:, 0]

    # for i, mass in enumerate(masses):
    #     # get the evolutionary track for all ages
    #     all_isochrones = []
    #     for age in all_ages_Gyr[all_ages_Gyr < 0.5]:
    #         all_isochrones.append(load_isochrone_generic(age * 1e3))
    #     myr_stack = np.vstack(all_isochrones).reshape(len(all_isochrones), 30, 11)

    #     mass_index = np.where(mass == mass_array)[0]
    #     mi = mass_index
    #     print(f"mass_index={mass_index}")

    #     teff_i = myr_stack[:, mi, 1]
    #     Mk_i = myr_stack[:, mi, -1]
    #     M_i = myr_stack[:, mi, 0]

    #     plt.plot(teff_i, Mk_i + distmod, "--", label=f"{mass} M$_{{\odot}}$")

    #     if M_i[0] <= 0.3:
    #         if M_i[0] in masslist:
    #             print("M = ", M_i[0])
    #             plt.plot(SpT_i, Mk_i + distmod, label=f"{M_i[0]} M$_{{\odot}}$")

    # ax.set_xlim(3, 13)
    # ax.invert_yaxis()
    # ax.invert_xaxis()
    ax.set_xlabel(r"$T_{\rm{eff}}$ (K)")
    ax.set_ylabel("Mass (M$_{\odot})$")
    plt.legend()

    return fig


def mock_lalchand_figure_v5():
    # Teff vs Mass!

    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)
    plt.grid(True, linestyle="--", lw=0.25)

    ages = [0.5, 2.0, 5.0, 10.0]

    for age, lw in zip(ages, lw_array):

        myr_x = load_isochrone_generic(age)

        mass = myr_x[:, 0]
        teff = myr_x[:, 1]
        Mk = myr_x[:, -1]
        M = myr_x[:, 0]

        u = M <= 0.8

        plt.plot(mass[u], teff[u], "k", lw=lw, label=f"{age} Myr")

    all_ages_Gyr = extract_age_array_Gyr(filepath)

    # masses = [0.2, 0.1, 0.05, 0.03, 0.01]
    # myr_0_5 = load_isochrone_generic(0.5)
    # mass_array = myr_0_5[:, 0]

    # for i, mass in enumerate(masses):
    #     # get the evolutionary track for all ages
    #     all_isochrones = []
    #     for age in all_ages_Gyr[all_ages_Gyr < 0.5]:
    #         all_isochrones.append(load_isochrone_generic(age * 1e3))
    #     myr_stack = np.vstack(all_isochrones).reshape(len(all_isochrones), 30, 11)

    #     mass_index = np.where(mass == mass_array)[0]
    #     mi = mass_index
    #     print(f"mass_index={mass_index}")

    #     teff_i = myr_stack[:, mi, 1]
    #     Mk_i = myr_stack[:, mi, -1]
    #     M_i = myr_stack[:, mi, 0]

    #     plt.plot(teff_i, Mk_i + distmod, "--", label=f"{mass} M$_{{\odot}}$")

    #     if M_i[0] <= 0.3:
    #         if M_i[0] in masslist:
    #             print("M = ", M_i[0])
    #             plt.plot(SpT_i, Mk_i + distmod, label=f"{M_i[0]} M$_{{\odot}}$")

    # ax.set_xlim(3, 13)
    # ax.invert_yaxis()
    # ax.invert_xaxis()
    ax.set_ylabel(r"$T_{\rm{eff}}$ (K)")
    ax.set_xlabel("Mass (M$_{\odot})$")
    plt.legend()

    return fig


def make_fig_3_HR_diagram():
    """
    Figure 3 is a combined HR diagram / observational summary of the objects in this study.

    """

    # This part pulls up the source information from our position-matching scheme.
    ngc_match = match_ngc()
    ic_match = match_ic()

    spread_dict = {
        "ngc": spreadsheet.load_wserv_v2(7),
        "ic": spreadsheet.load_wserv_v2(8),
    }
    match_dict = {"ngc": ngc_match, "ic": ic_match}
    fullname_dict = {"ngc": "NGC 1333", "ic": "IC 348"}
    names = ["ic", "ngc"]

    # This part sets up the figure.

    plt.rcParams["font.size"] = 6

    fig = plt.figure(figsize=(8, 8), dpi=300)

    # There will be eight subplots. Let's prototype.

    # cmd_axes = [fig.add_subplot(2, 2, 1), fig.add_subplot(2, 2, 2)]
    # hist_axes = [fig.add_subplot(6, 2, 6 + 1), fig.add_subplot(6, 2, 6 + 2)]
    # hr_axes = [fig.add_subplot(6, 2, 8 + 1), fig.add_subplot(6, 2, 8 + 2)]
    # sigma_axes = [fig.add_subplot(6, 2, 10 + 1), fig.add_subplot(6, 2, 10 + 2)]

    # CMD
    cmd1 = fig.add_subplot(2, 2, 1)
    cmd2 = fig.add_subplot(2, 2, 2, sharex=cmd1, sharey=cmd1)
    cmd_axes = [cmd1, cmd2]

    # Hist
    hist1 = fig.add_subplot(6, 2, 7)
    hist2 = fig.add_subplot(6, 2, 8, sharex=hist1, sharey=hist1)
    hist_axes = [hist1, hist2]

    # HR
    hr1 = fig.add_subplot(6, 2, 9)
    hr2 = fig.add_subplot(6, 2, 10, sharex=hr1, sharey=hr1)
    hr_axes = [hr1, hr2]

    # Sigma
    sigma1 = fig.add_subplot(6, 2, 11)
    sigma2 = fig.add_subplot(6, 2, 12, sharex=sigma1, sharey=sigma1)
    sigma_axes = [sigma1, sigma2]

    # (This'll be a loop, eventually)
    for i in range(len(names)):
        name = names[i]
        cmd_ax = cmd_axes[i]
        hist_ax = hist_axes[i]
        hr_ax = hr_axes[i]
        sigma_ax = sigma_axes[i]

        fullname = fullname_dict[name]
        match = match_dict[name]
        spread = spread_dict[name]
        ir_exc = match.approved["IRexc"] == "yes"

        ################
        # Draw the CMD #
        ################

        # 1. the data

        cmd_ax.plot(
            match.approved["median_HMKPNT"][~ir_exc],
            match.approved["median_KAPERMAG3"][~ir_exc],
            "ks",
            ms=1,
            label="no IR exc",
        )
        cmd_ax.plot(
            match.approved["median_HMKPNT"][ir_exc],
            match.approved["median_KAPERMAG3"][ir_exc],
            "rD",
            ms=1,
            label="with IR exc",
        )
        # cmd_ax.plot(
        #     spread["median"]["HMKPNT"],
        #     spread["median"]["KAPERMAG3"],
        #     "k,",
        #     alpha=0.1,
        #     zorder=-1,
        # )
        cmd_ax.set_title(f"{fullname}")
        cmd_ax.set_xlabel("$H-K$ color", labelpad=-0.04)
        cmd_ax.set_ylabel("$K$ mag")
        # cmd_ax.set_xlim(0, 4.05)
        cmd_ax.set_xlim(0, 2.5)
        cmd_ax.set_ylim(18.5, 9)

        if i == 0:
            dot_legend = cmd_ax.legend(bbox_to_anchor=(1.0, 0.85), loc="upper right")
            cmd_ax.add_artist(dot_legend)

        # 2. the theory

        lws = [0.2, 0.75]
        ages = [1, 8]

        isochrone_lines = []
        for age, lw in zip(ages, lws):

            myr_x = load_isochrone_generic(age)

            M_k = myr_x[:, 10]
            color = myr_x[:, 9] - myr_x[:, 10]

            isochrone_lines.append(
                cmd_ax.plot(color, M_k + distmod, "k-", lw=lw, label=f"{age} Myr")[0]
            )

            # But let's also do three masses... 0.5, 0.08, 0.01

        #
        masses = [0.01, 0.08, 0.2][::-1]
        myr_0_5 = load_isochrone_generic(0.5)
        mass_array = myr_0_5[:, 0]
        lw_array = [0.5, 1, 1.5][::-1]

        all_ages_Gyr = extract_age_array_Gyr(filepath)

        for j, (mass, lw) in enumerate(zip(masses, lw_array)):
            # get the evolutionary track for all ages
            all_isochrones = []
            for age in all_ages_Gyr[all_ages_Gyr < 0.5]:
                all_isochrones.append(load_isochrone_generic(age * 1e3))
            myr_stack = np.vstack(all_isochrones).reshape(len(all_isochrones), 30, 11)

            mass_index = np.where(mass == mass_array)[0]
            mi = mass_index

            Mk_i = myr_stack[:, mi, -1]
            Mh_i = myr_stack[:, mi, -2]

            isochrone_lines.append(
                cmd_ax.plot(
                    Mh_i - Mk_i,
                    Mk_i + distmod,
                    "C0--",
                    lw=lw,
                    label=f"{mass} M$_{{\odot}}$",
                )[0]
            )

        # plt.gca().invert_yaxis()
        if i == 0:

            age_legend = cmd_ax.legend(
                handles=isochrone_lines,
                labels=[x.get_label() for x in isochrone_lines],
                title="$d=300$ pc",
                loc="lower right",
                bbox_to_anchor=(-0.05, 0.0),
            )
            cmd_ax.add_artist(age_legend)

        ######################
        # Draw the histogram #
        ######################

        hist_ax.hist(
            match.approved["SpT"],
            range=[0, 15],
            bins=np.arange(0, 15, 0.5),
            edgecolor="k",
            linewidth=0.5,
            facecolor="None",
            histtype="stepfilled",
            label="all targets: " f"$n$={len(match.approved)}",
        )
        hist_ax.hist(
            match.statistical["SpT"],
            edgecolor="k",
            facecolor="0.5",
            # hatch="..",
            range=[0, 15],
            bins=np.arange(0, 15, 0.5),
            histtype="stepfilled",
            label="high-quality photometry \n($\sigma_{K} \\leq 0.05$ mag): "
            f"$n$={len(match.statistical)}",
        )

        hist_ax.set_xlim(0, 14)
        hist_ax.set_ylabel("Number of sources")
        hist_ax.set_xlabel("Spectral Type")
        hist_ax.legend()

        spt_array = np.array([get_SpT_from_num(int(x)) for x in hist_ax.get_xticks()])
        hist_ax.xaxis.set_tick_params(labelbottom=True)
        hist_ax.set_xticklabels(spt_array)

        hist_ax.grid(True, axis="x", ls=":")

        #######################
        # Draw the HR diagram #
        #######################

        hr_ax.plot(
            match.approved["SpT"], match.approved["median_KAPERMAG3"], "k.", ms=2
        )
        hr_ax.set_xlim(0, 14)
        hr_ax.set_ylabel("median $K$ mag")
        hr_ax.set_xlabel("Spectral Type")
        hr_ax.grid(True, axis="x", ls=":")

        spt_array = np.array([get_SpT_from_num(int(x)) for x in hr_ax.get_xticks()])
        hr_ax.xaxis.set_tick_params(labelbottom=True)
        hr_ax.set_xticklabels(spt_array)

        if i == 0:
            hr_ax.invert_yaxis()

        sigma_ax.plot(
            match.approved["SpT"], match.approved["median_KAPERMAG3ERR"], "k.", ms=2
        )
        sigma_ax.set_xlim(0, 14)
        sigma_ax.semilogy()
        sigma_ax.set_xlabel("Spectral Type")
        sigma_ax.set_ylabel("$\sigma_{K}$ (mag)")

        # ax.set_xticklabels(spt_array)
        sigma_ax.set_xticklabels(spt_array)
        sigma_ax.grid(True, axis="x", ls=":")
        sigma_ax.axhline(0.05, lw=0.5, ls="--")

    plt.show()
    fig.savefig(
        os.path.join(figure_export_path, "Figure_3_CMD_HR_diagram.pdf"),
        bbox_inches="tight",
    )

    return fig


if __name__ == "__main__":
    # mock_lalchand_figure()
    fig = make_fig_3_HR_diagram()
