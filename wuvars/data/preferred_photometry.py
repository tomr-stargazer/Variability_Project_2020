"""

A module that loads the preferred photometry for easy import.

"""

from wuvars.data.photometry import group_wserv_v2, load_wserv_v2

photometry_wserv7 = group_wserv_v2(load_wserv_v2(7, suffix="_outlier_cleaned_152"))

photometry_wserv8 = group_wserv_v2(load_wserv_v2(8))
