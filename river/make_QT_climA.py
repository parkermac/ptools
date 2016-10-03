#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 09:13:43 2016

@author: PM5

Program to gather create climatological records for rivers.  Does both
flow and temperature.

** Analytical version **
"""

#%% imports
import os
import sys

alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)

from importlib import reload
import Lfun
reload(Lfun)
Ldir = Lfun.Lstart(gridname='test')

import pandas as pd
import numpy as np

#%% create directory for output, if needed
dir0 = Ldir['data'] + 'rivers/'
Lfun.make_dir(dir0, clean=False)

out_dir = dir0 + 'Data_clim/'
Lfun.make_dir(out_dir, clean=False)

out_dir_T = dir0 + 'Data_T_clim/'
Lfun.make_dir(out_dir_T, clean=False)


#%% create the climatologies

rn = 'creek1'

# flow
qtc = pd.Series(dict(zip(range(1,367), 500. * np.ones(366))))
qtc.to_csv(out_dir + rn + '.csv')

# temperature
ttc = pd.Series(dict(zip(range(1,367), 10. * np.ones(366))))
ttc.to_csv(out_dir_T + rn + '.csv')

