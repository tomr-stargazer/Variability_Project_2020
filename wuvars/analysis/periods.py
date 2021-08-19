"""
Some common utilities to help us with period-finding.

References the `Following the recommendations of VanderPlas` notebook.

"""

import numpy as np

f_min = 0  # no longest period
f_max = 24  # 24 day^-1 (i.e., shortest period is 1 hr), due to astrophysical prior on shortest rotation rate observed for young brown dwarfs (Vos et al. 2020)
N_eval = int(2.4e5)  # (240,000 evaluations, which oversamples each 'peak' by about a factor 20, very similar to the recommendation in Section 7.1 of vanderPlas 2018)

freq = np.linspace(f_min, f_max, N_eval)
