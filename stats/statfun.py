# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 09:55:33 2016

@author: PM5

Functions for statistics.

"""

import numpy as np

def cross_correlation(x,y):

    """
    Returns the cross-correlation coefficient for a pair of time series, x
    and y.

    It is defined as in Oke et al. Part I (2002) JPO, vol 32, p. 1379,
    Eqn. (C1)
    """

    # first flatten the records
    x = x.flatten()
    y = y.flatten()
    # then remove all nans
    mask = (x!=np.nan) & (y!=np.nan)
    x = x[mask]
    y = y[mask]

    if len(x)>0 and len(y)>0:
        cc = np.mean( (x-np.mean(x)) * (y-np.mean(y)) ) / (np.std(x) * np.std(y))
    else:
        cc = np.nan

    return cc
