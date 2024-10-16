"""
An importable set of magnitudes of brown dwarf limit stars.

These are from an email Johanna Vos sent me on July 21, 2020. She wrote:
"This text file contains absolute and apparent mags (J,H,K) for mass<0.08Msolar for 
each star-forming region. I used the Allard BT-Settl models to calculate these. 
I didn't do any sort of error calculation but can implement an MCMC type 
error propagation if necessary."

"""

absolute_BD_mags_jhk = {}
apparent_BD_mags_jhk = {}

# ------------------------------------------------------
# Cygnus OB7

# Absolute 2MASS mags:
absolute_BD_mags_jhk[1] = (7.8654, 7.3, 7.0216)
# Apparent 2MASS mags:
apparent_BD_mags_jhk[1] = (17.38084993495972, 16.81544993495972, 16.537049934959718)
# ------------------------------------------------------

# ------------------------------------------------------
# Orion Nebula Cluster

# Absolute 2MASS mags:
absolute_BD_mags_jhk[5] = (6.063, 5.5, 5.234)
# Apparent 2MASS mags:
apparent_BD_mags_jhk[5] = (14.148001705604496, 13.585001705604496, 13.319001705604496)
# ------------------------------------------------------

# ------------------------------------------------------
# NGC 1333

# Absolute 2MASS mags:
absolute_BD_mags_jhk[7] = (6.0145, 5.4535, 5.186)
# Apparent 2MASS mags:
apparent_BD_mags_jhk[7] = (13.348838101770548, 12.787838101770548, 12.520338101770548)
# ------------------------------------------------------

# ------------------------------------------------------
# IC348

# Absolute 2MASS mags:
absolute_BD_mags_jhk[8] = (7.155, 6.59, 6.317)
# Apparent 2MASS mags:
apparent_BD_mags_jhk[8] = (14.687525162024361, 14.12252516202436, 13.84952516202436)
# ------------------------------------------------------

# ------------------------------------------------------
# Mon R2

# Absolute 2MASS mags:
absolute_BD_mags_jhk[11] = (6.952, 6.386, 6.114)
# Apparent 2MASS mags:
apparent_BD_mags_jhk[11] = (16.706257294442732, 16.140257294442733, 15.868257294442731)
# ------------------------------------------------------


# Johanna's `BD_magnitudes.txt` from July 21, 2020:

# ------------------------------------------------------
# Cygnus OB7

# Absolute 2MASS mags:
# 7.8654 7.3 7.0216
# Apparent 2MASS mags:
# 17.38084993495972 16.81544993495972 16.537049934959718
# ------------------------------------------------------

# ------------------------------------------------------
# Orion Nebula Cluster

# Absolute 2MASS mags:
# 6.063 5.5 5.234
# Apparent 2MASS mags:
# 14.148001705604496 13.585001705604496 13.319001705604496
# ------------------------------------------------------

# ------------------------------------------------------
# NGC 1333

# Absolute 2MASS mags:
# 6.0145 5.4535 5.186
# Apparent 2MASS mags:
# 13.348838101770548 12.787838101770548 12.520338101770548
# ------------------------------------------------------

# ------------------------------------------------------
# IC348

# Absolute 2MASS mags:
# 7.155 6.59 6.317
# Apparent 2MASS mags:
# 14.687525162024361 14.12252516202436 13.84952516202436
# ------------------------------------------------------

# ------------------------------------------------------
# Mon R2

# Absolute 2MASS mags:
# 6.952 6.386 6.114
# Apparent 2MASS mags:
# 16.706257294442732 16.140257294442733 15.868257294442731
# ------------------------------------------------------
