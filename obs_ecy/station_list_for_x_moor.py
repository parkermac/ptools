"""
This makes a list of all the Ecology and EC CTD stations in a format
suitable for feeding to x_moor.
"""

import pandas as pd
import numpy as np
import pickle

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
Ldir = Lfun.Lstart()

# load processed station info
dir0 = Ldir['parent'] + 'ptools_data/ecology/'
sta_df0 = pd.read_pickle(dir0 + 'sta_df.p')
dir1 = Ldir['parent'] + 'ptools_data/canada/'
sta_df1 = pd.read_pickle(dir1 + 'sta_df.p')
sta_df = pd.concat((sta_df0, sta_df1), sort=False)
sta_df = sta_df[['Latitude', 'Longitude']]

# print to screen for copying to x_moor/moor_lists.py
for sta in sta_df.index:
    x = sta_df.loc[sta,'Longitude']
    y = sta_df.loc[sta,'Latitude']
    print("\'%s\': (%0.6f, %0.6f)," % (sta, x, y))
