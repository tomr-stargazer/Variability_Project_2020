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
from wuvars.analysis.spectral_type_to_number import (get_num_from_SpT,
                                                     get_SpT_from_num)
from wuvars.analysis.spectral_type_to_temperature import get_SpT_from_Teff_HH14

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


lw_array = [0.5, 0.75, 1.0, 1.25]

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

    fig, ax = plt.subplots(1, dpi=150)

    ages = [0.5, 2.0, 5.0, 10.0]

    for age, lw in zip(ages, lw_array):

        myr_x = load_isochrone_generic(age)

        mass = myr_x[:, 0]
        teff = myr_x[:, 1]
        Mk = myr_x[:, -1]
        M = myr_x[:, 0]

        u = M <= 0.3

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
    ax.set_ylabel("Mass (M$_{\odot}$")
    plt.legend()

    return fig

if __name__ == "__main__":
    mock_lalchand_figure()
