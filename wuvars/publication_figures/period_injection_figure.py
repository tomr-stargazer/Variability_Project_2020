"""
I'm generating a period-injection figure.

"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from make_figures import figure_export_path

save_path = (
    "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes"
)

periods_denser_raw = np.logspace(np.log10(1 / 24.5), np.log10(24.5), 60)
extra_periods = (1 / 4, 1 / 3, 1 / 2, 1, 2, 3, 4)
periods_denser = np.sort(np.append(periods_denser_raw, np.array(extra_periods)))

amplitudes_denser = [
    0.001,
    0.003,
    0.005,
    0.0075,
    0.0085,
    0.01,
    0.0125,
    0.015,
    0.0175,
    0.02,
    0.025,
    0.03,
]

output_reconstructed = []

for i in range(4):
    a = np.load(os.path.join(save_path, f"output_denser_{i}_20210823.npy"))
    output_reconstructed.append(a)

output_denser = output_reconstructed


def make_period_injection_fig():

    plt.style.use("default")

    mpl.rcParams["xtick.direction"] = "out"
    mpl.rcParams["ytick.direction"] = "out"
    mpl.rcParams["xtick.top"] = False
    mpl.rcParams["ytick.right"] = False

    fig, ax = plt.subplots(1, figsize=(5.5, 4), dpi=150)
    plt.contourf(
        periods_denser,
        amplitudes_denser,
        100 * np.array(output_denser[0]).T,
        levels=100 * np.linspace(0, 1, 11),
        cmap="RdBu",
    )
    plt.ylim(0, 0.03)
    plt.xlabel("Injected period (days)")
    plt.ylabel("Injected amplitude (mag)")
    # plt.title(r"Period recovery accuracy ($K < 15$, n=40)")
    plt.semilogx()
    # plt.xlim(None,2.6)
    cbar = plt.colorbar(pad=0.02)
    cbar.set_label("Periods correctly recovered (%)", labelpad=-2.5)

    for per in extra_periods:
        plt.axvline(per, lw=0.75, ls="--", color="w", dashes=(5, 10))

    plt.plot(
        [0, 0],
        [0, 0],
        lw=0.75,
        ls="--",
        color="w",
        dashes=(5, 10),
        label="common\naliases",
    )
    plt.legend(loc="upper right", fontsize=8, framealpha=0.5)

    # fig.savefig("period_recovery_grid_WSERV8.png")
    # fig.savefig("period_recovery_grid_WSERV8.pdf")
    fig.savefig(os.path.join(figure_export_path, "Figure_5_Period_Recovery.pdf"))

    return fig


if __name__ == "__main__":
    make_period_injection_fig()
