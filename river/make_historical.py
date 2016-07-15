"""
Program to gather historical records for rivers.
"""

#%% imports
import matplotlib.pyplot as plt
import os
import sys

alp = os.path.abspath('../../alpha')
if alp not in sys.path:
    sys.path.append(alp)


from importlib import reload
import Lfun
reload(Lfun)
Ldir = Lfun.Lstart(gridname='test')

rivp = os.path.abspath(Ldir['LO'] + 'forcing/riv1/')
if rivp not in sys.path:
    sys.path.append(rivp)

import river_class
reload(river_class)

from datetime import datetime
import pandas as pd
import numpy as np

#%% Load a dataframe with info for rivers to get

ri_fn = Ldir['parent'] + 'ptools_output/river/pnw_all_2016_07/river_info.csv'

df = pd.read_csv(ri_fn, index_col='rname')

# and create directory for output, if needed
out_dir0 = Ldir['data'] + 'rivers/'
out_dir = out_dir0 + 'Data_historical/'
Lfun.make_dir(out_dir0, clean=False)
Lfun.make_dir(out_dir, clean=False)

#%% get river info

dt0 = datetime(2015,1,1)
dt1 = datetime(2015,12,31)
days = (dt0, dt1)

qt_dict = dict()

for rn in df.index:

    if rn == 'fraser':

        rs = df.ix[rn] # a series with info for this river
        riv = river_class.River(rn, rs, days)

        if not pd.isnull(rs.usgs):
            riv.get_usgs_data(days)
        elif not pd.isnull(rs.ec) and dt1 >= datetime(2015,1,1):
            riv.get_ec_data(days)
        elif not pd.isnull(rs.ec) and dt1 <= datetime(2014,12,31):
            # always gives a year of data, for dt0.year
            riv.get_ec_data_historical(days)

        riv.print_info()
        sys.stdout.flush()

        if not riv.qt.empty:
            qt_dict[rn] = riv.qt

#%% save output

for rn in qt_dict.keys():

    qt = qt_dict[rn]

    qt.to_pickle(out_dir + rn + '.p')

#%% plotting

if True:

    plt.close()

    NP = len(qt_dict)

    NR = np.maximum(1, np.ceil(np.sqrt(NP)).astype(int))
    NC = np.ceil(np.sqrt(NP)).astype(int)

    fig, axes = plt.subplots(nrows=NR, ncols=NC, figsize=(17,9), squeeze=False)

    cc = 0
    for rn in qt_dict.keys():
        ir = int(np.floor(cc/NC))
        ic = int(cc - NC*ir)
        ax = axes[ir, ic]
        qt_dict[rn].plot(ax=ax, title=rn, style='-k')
        cc += 1

    plt.show()


