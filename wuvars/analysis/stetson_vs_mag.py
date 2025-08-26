"""
This is a script which generates the Stetson threshold vs magnitude plots.

"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wuvars.analysis.variability_selection as vs
from wuvars.analysis.moving_average_variability_threshold import (
    moving_average_percentile, moving_average_percentile_monotonic)
from wuvars.analysis.stetson_vs_mag_selection_fns import (ic_selection_fn,
                                                          ngc_selection_fn,
                                                          onc_selection_fn)
from wuvars.data import spreadsheet

spread = spreadsheet.load_v2()

selection_fns = {}
selection_fns[5] = onc_selection_fn
selection_fns[7] = ngc_selection_fn
selection_fns[8] = ic_selection_fn

wserv_ids = [1, 5, 7, 8, 11]
n_min_list = [60, 35, 80, 55, 65]
n_max_list = [100, 90, 160, 80, 100]

min_Stetson_list = [2, 3, 1.6, 1.6, 4]
SFR_names = ["Cyg OB7", "Orion Nebula Cluster", "NGC 1333", "IC 348", "Mon R2"]
SFR_dict = {x: y for x, y in zip(wserv_ids, SFR_names)}


if __name__ == "__main__":

    dpi = 150

    for wserv, n_min, n_max, S in list(
        zip(wserv_ids, n_min_list, n_max_list, min_Stetson_list)
    )[1:-1]:

        ds = getattr(spread, f"wserv{wserv}")

        q0 = vs.sq0(ds, n_min, n_max)
        q1 = vs.sq1(ds, n_min, n_max)
        q2 = vs.sq2(ds, n_min, n_max)

        v0 = vs.sq0_variables(ds, n_min, n_max, Stetson_cutoff=S)
        v1 = vs.sq1_variables(ds, n_min, n_max, Stetson_cutoff=S)
        v2 = vs.sq2_variables(ds, n_min, n_max, Stetson_cutoff=S)

        selection_region = selection_fns[wserv](ds)

        ### first figure: simply a map of what region we are using as a 'control'.

        fig, ax = plt.subplots(ncols=1, figsize=(6, 6), dpi=dpi)

        ax.plot(
            np.degrees(ds[q0 & ~selection_region]["mean"]["RA"]),
            np.degrees(ds[q0 & ~selection_region]["mean"]["DEC"]),
            "k,",
            alpha=0.1,
        )
        ax.invert_xaxis()

        ax.plot(
            np.degrees(ds[v1 & ~selection_region]["mean"]["RA"]),
            np.degrees(ds[v1 & ~selection_region]["mean"]["DEC"]),
            "r.",
        )

        fig.suptitle(f"WSERV{wserv}: {SFR_dict[wserv]}")

        ax.set_title("Variable stars (red)")
        ax.set_xlabel("R.A. (deg)")
        ax.set_ylabel("Decl. (deg)")

        fig2, ax2 = plt.subplots(ncols=1, figsize=(6, 6), dpi=dpi)

        ax2.plot(
            ds[q2 & ~selection_region]["median"]["KAPERMAG3"],
            ds[q2 & ~selection_region]["variability"]["Stetson_JHK"],
            "k.",
            ms=2,
        )

        xs = ds[q2 & ~selection_region]["median"]["KAPERMAG3"]
        ys = ds[q2 & ~selection_region]["variability"]["Stetson_JHK"]
        x_grid, result_grid_50 = moving_average_percentile(
            xs, ys, x_range=(11, 18), percentile=50
        )
        x_grid, result_grid_95 = moving_average_percentile_monotonic(
            xs, ys, x_range=(11, 18), percentile=95
        )
        x_grid, result_grid_997 = moving_average_percentile_monotonic(
            xs, ys, x_range=(11, 18), percentile=99.7
        )

        save_dir = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "stetson_v_mag"
        )
        np.save(
            os.path.join(save_dir, f"WSERV{wserv}_result_grid_997"),
            np.vstack([x_grid, result_grid_997]),
        )
        np.save(
            os.path.join(save_dir, f"WSERV{wserv}_result_grid_95"),
            np.vstack([x_grid, result_grid_95]),
        )
        np.save(
            os.path.join(save_dir, f"WSERV{wserv}_result_grid_50"),
            np.vstack([x_grid, result_grid_50]),
        )



        ax2.plot(
            ds[q2 & v2 & ~selection_region]["median"]["KAPERMAG3"],
            ds[q2 & v2 & ~selection_region]["variability"]["Stetson_JHK"],
            "r.",
        )
        ax2.plot(x_grid, result_grid_50, "C0-.", lw=2, label=r"median")
        ax2.plot(x_grid, result_grid_95, "b--", lw=2, label=r"2$\sigma$")
        ax2.plot(x_grid, result_grid_997, "r-", lw=2, label=r"3$\sigma$")

        ax2.semilogy()
        ax2.axhline(S, color="k", lw=1, alpha=0.5)
        ax2.set_ylim(1e-2, None)

        ax2.set_xlim(None, max(ds[q2 & ~selection_region]["median"]["KAPERMAG3"]))
        ax2.set_ylabel("Stetson index (JHK)")
        ax2.set_xlabel("median $K$ magnitude")

        ax2.set_title(
            f"WSERV{wserv} ({SFR_dict[wserv]}) Stetson variability index vs. $K$ mag"
        )

        ax2.legend()
        fig.savefig(
            os.path.join(save_dir, f"WSERV{wserv}_map_selection.png"),
        )
        fig2.savefig(
            os.path.join(save_dir, f"WSERV{wserv}_Stetson_v_mag.png"),
        )
        plt.show()
