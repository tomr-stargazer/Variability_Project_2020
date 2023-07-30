"""
In this script I am coding up maps of the regions we have studied. 

I want to set this up for publication.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from wuvars.data import spreadsheet, photometry, quality_classes
from wuvars.analysis.bd_matching_ngc1333 import lowmass_ngc_joint_matches

lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"
ngc_spread = spreadsheet.load_wserv_v2(7)
ngc_q = quality_classes.load_q(7)
inspect_ngc = pd.read_csv(
    os.path.join(lc_dir, "inspection_ngc.csv"), skipinitialspace=True
)


# Step one: reproduce literally what has been done.

# I am using the following file:
# ./wuvars/analysis/prototypes/Lightcurves of all matches - first pass for inspection.ipynb

plt.figure(figsize=(8, 8), dpi=150)

rejected_sources_ngc = inspect_ngc["SOURCEID"][inspect_ngc["exclude?"] == "yes"]
# inspect_onc['SOURCEID']

rejected_sources_ngc_ind = np.in1d(
    lowmass_ngc_joint_matches["SOURCEID"], rejected_sources_ngc
)

plt.plot(
    np.degrees(ngc_spread["median"]["RA"]),
    np.degrees(ngc_spread["median"]["DEC"]),
    "k,",
    alpha=0.5,
    rasterized=True,
)

plt.plot(
    np.degrees(lowmass_ngc_joint_matches["RA"][~rejected_sources_ngc_ind]),
    np.degrees(lowmass_ngc_joint_matches["DEC"][~rejected_sources_ngc_ind]),
    "k.",
    label="Luhman LMS+BDs with usable lightcurve data (n=103)",
)
plt.plot(
    np.degrees(lowmass_ngc_joint_matches["RA"][rejected_sources_ngc_ind]),
    np.degrees(lowmass_ngc_joint_matches["DEC"][rejected_sources_ngc_ind]),
    "rx",
    ms=6,
    mew=1.7,
    label="Luhman LMS+BDs with unusable lightcurve data (n=2)",
)
plt.gca().invert_xaxis()
plt.xlabel("RA (deg)")
plt.ylabel("Dec (deg)")

plt.xlim(
    np.degrees(ngc_spread["median"]["RA"][ngc_q.q2]).max(),
    np.degrees(ngc_spread["median"]["RA"][ngc_q.q2]).min(),
)
plt.ylim(
    np.degrees(ngc_spread["median"]["DEC"]).min(),
    np.degrees(ngc_spread["median"]["DEC"]).max(),
)

plt.legend(fontsize=14, loc="lower left")
plt.title("Accepted and rejected source matches in NGC 1333")

plt.savefig("test_NGC1333_map.png")
plt.savefig("test_NGC1333_map.pdf")

# === IC 348 ===

from wuvars.analysis.bd_matching_ic348 import lowmass_ic_joint_matches

lc_dir = "/Users/tsrice/Documents/Variability_Project_2020/wuvars/analysis/prototypes/BD_lcs_v3"
inspect_ic = pd.read_csv(
    os.path.join(lc_dir, "inspection_ic.csv"), skipinitialspace=True
)

ic_spread = spreadsheet.load_wserv_v2(8)
ic_q = quality_classes.load_q(8)

plt.figure(figsize=(8, 8), dpi=150)

rejected_sources_ic = inspect_ic["SOURCEID"][inspect_ic["exclude?"] == "yes"]
# inspect_onc['SOURCEID']

rejected_sources_ic_ind = np.in1d(
    lowmass_ic_joint_matches["SOURCEID"], rejected_sources_ic
)

plt.plot(
    np.degrees(ic_spread["median"]["RA"]),
    np.degrees(ic_spread["median"]["DEC"]),
    "k,",
    alpha=0.5,
    rasterized=True,
)

plt.plot(
    np.degrees(lowmass_ic_joint_matches["RA"][~rejected_sources_ic_ind]),
    np.degrees(lowmass_ic_joint_matches["DEC"][~rejected_sources_ic_ind]),
    "k.",
    label="Luhman LMS+BDs with usable lightcurve data (n=225)",
)
plt.plot(
    np.degrees(lowmass_ic_joint_matches["RA"][rejected_sources_ic_ind]),
    np.degrees(lowmass_ic_joint_matches["DEC"][rejected_sources_ic_ind]),
    "rx",
    ms=6,
    mew=1.7,
    label="Luhman LMS+BDs with unusable lightcurve data (n=9)",
)
plt.gca().invert_xaxis()
plt.xlabel("RA (deg)")
plt.ylabel("Dec (deg)")

plt.xlim(
    np.degrees(ic_spread["median"]["RA"][ic_q.q2]).max(),
    np.degrees(ic_spread["median"]["RA"][ic_q.q2]).min(),
)
plt.ylim(
    np.degrees(ic_spread["median"]["DEC"]).min(),
    np.degrees(ic_spread["median"]["DEC"]).max(),
)

plt.legend(fontsize=14, loc="lower left")
plt.title("Accepted and rejected source matches in IC 348")

plt.savefig("test_IC348_map.png")
plt.savefig("test_IC348_map.pdf")
