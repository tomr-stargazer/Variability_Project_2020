"""
I plan to put some light curve plotting code here.

"""

import numpy as np
import matplotlib.pyplot as plt


def simple_lc(dat):

    # set up data

    date = dat['MEANMNJDOBS']

    j = dat['JAPERMAG3']
    h = dat['HAPERMAG3']
    k = dat['KAPERMAG3']

    j_e = dat['JAPERMAG3ERR']
    h_e = dat['HAPERMAG3ERR']
    k_e = dat['KAPERMAG3ERR']


    # set up plot

    fig = plt.figure(figsize = (10, 6), dpi=80,
                     facecolor='w', edgecolor='k')

    bottom = 0.1
    height = .25
    left = 0.075
    width = 0.5

    ax_k = fig.add_axes( (left, bottom, width, height) )
    ax_h = fig.add_axes( (left, bottom+.3, width, height), sharex=ax_k )
    ax_j = fig.add_axes( (left, bottom+.6, width, height), sharex=ax_k )
    
    ax_jhk = fig.add_axes( (.65, bottom, .3, .375) )
    ax_khk = fig.add_axes( (.65, bottom+.475, .3, .375) )

    ax_j.errorbar(date, j, yerr=j_e, fmt='bo', ecolor='k', ms=2)
