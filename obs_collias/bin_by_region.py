#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Next steps in processing data:

(1) collect everything in regions

Regions (the first digit is stored ar Num0 in Bottles):
        100        Strait of Juan de Fuca
        200        Admiralty Inlet
        300        Puget Sound Basin
        400        Southern Puget Sound
        500        Hood Canal
        600        Whidbey Basin
        700        North Sound
        800        San Juan Island Passages

(2) and interpolate to selected depths, typically -30, -10, and 0 m
to match the bottle depths from the modern Ecology data.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun

pth = os.path.abspath(Ldir['LO'] + 'plotting')
if pth not in sys.path:
    sys.path.append(pth)
import pfun

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

dir0 = Ldir['parent'] + 'ptools_data/collias/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

dir1 = Ldir['parent'] + 'ptools_output/collias/'
Lfun.make_dir(dir1)

testing = False
if testing:
    year_list = [1952]
else:
    year_list = range(1932, 1976) # data go from 1932 to 1975

# variables to process
col_list = ['Temp. (deg C)', 'Salinity',
            'DO (mg L-1)', 'NO3 (mg L-1)', 'NO2 (mg L-1)', 'SiOH4 (mg L-1)']

for region in range(1,9):
    print('Binning by region ' + str(region))
    a = pd.DataFrame()
    for year in year_list:
        print(' - year = ' + str(year))
        b = pd.read_pickle(dir0 + 'Bottles_' + str(year) + '.p')
    
        if len(b) > 0:
            b = b[b['Num0']==str(region)] # select region
    
            # now we want to go through all individual casts
            sta_list = b['Station'].values
            for sta in set(sta_list):
                bb = b[b['Station']==sta]
                date_list = bb['Date'].values
                for date in set(date_list):
                    c = bb[bb['Date']==date]
                    c = c.set_index('Z (m)')
                    c = c.reindex(columns=col_list)
                    # now we have a litte DataFrame that is a single cast and
                    # so this is a good place to interpolate to just selected depths
                    zvec = np.array(c.index.values)
                    z = np.array([-30, -10, 0])
                    if len(zvec) > 1:
                        ii0, ii1, ffr = zfun.get_interpolant(z,zvec, extrap_nan=True)
                        x = pd.DataFrame(columns=col_list) #  a temporary storage place
                        # NOTE the brackets around [i0 or i1] keep the result as a DataFrame even
                        # though we are only pulling out a single row.
                        iz = 0
                        for fr in ffr:
                            i0 = ii0[iz]; i1 = ii1[iz]; zz = z[iz]
                            d = (1-fr)*c.iloc[[i0],:].values + fr*c.iloc[[i1],:].values
                            x.loc[zz,:] = d#dict(zip(col_list,d))
                            iz+=1
                        x.loc[:,'Date'] = date
                        a = pd.concat((a,x))

    if len(a) > 0:
        a.loc[:,'Z (m)'] = a.index.values
        a = a.set_index('Date')
        a.to_pickle(dir1 + 'region_' + str(region) + '.p')

