# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 08:43:12 2016

@author: PM5

Test of time conversions.

"""

from datetime import datetime, timedelta
import numpy as np

def dt64_to_dt(dt64):
    # convert numpy datetime64 to datetime
    dt = datetime.utcfromtimestamp(dt64.astype('datetime64[ns]').tolist()/1e9)
    return dt

def dn_to_dt(dn):
    # convert MATLAB datenum to datetime
    # there is a tiny roundoff error involved, amounts to 5e-6 sec
    dt = ( datetime.fromordinal(int(dn))
        + timedelta(days = np.mod(dn,1))
        - timedelta(days = 366) )
    return dt

def dt_to_dn(dt):
    # convert datetime to MATLAB datenum
    # there is a tiny roundoff error involved, amounts to 1e-5 sec
    dn = ( datetime.toordinal(dt)
        + dt.hour/24 + dt.minute/(24*60) + dt.second/(86400)
        + dt.microsecond/(86400*1e6)
        + 366 )
    return dn

# generate original data

dt = datetime(2015, 7, 3, 11, 47, 2)

# dn was generated in MATLAB from the commands
# dn = datenum(2015, 7, 3, 11, 47, 2); format long
dn = 7.361484909953703e+05

# this converts from datetime to datetime64, so we don't need
# to make a special function to do it
dt64 = np.datetime64(dt)

# make conversions

dt1 = dt64_to_dt(dt64)

dt2 = dn_to_dt(dn)

dn1 = dt_to_dn(dt)


