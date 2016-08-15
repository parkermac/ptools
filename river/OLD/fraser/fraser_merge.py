# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 14:53:27 2016

@author: PM5

This merges the Fraser Series, which have some overlap.
"""

import os
import sys
pth = os.path.abspath('../../alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

import pandas as pd

indir = Ldir['data'] + 'rivers/data_processed/'
a = pd.read_pickle(indir + 'fraser_flow_historical.p')
b = pd.read_pickle(indir + 'fraser_flow_historical_2.p')
c = pd.read_pickle(indir + 'fraser_flow.p')

d = pd.concat([a,b,c])
dd = dd = d[~d.index.duplicated(keep='last')]
dd = dd.sort_index()

dd.plot()

dd.to_pickle(indir + 'fraser_flow_merged.p')

