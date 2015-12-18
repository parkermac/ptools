"""
Program to create long river records and analyze them.
"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
from importlib import reload
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart('cascadia1','base')

pth = os.path.abspath('../../LiveOcean/forcing/riv2')
if pth not in sys.path: sys.path.append(pth)
import river_class; reload(river_class)

import pandas as pd

outdir = Ldir['data'] + 'rivers/data_processed/'

if True:
    # get the list of rivers that we need for a run
    rdf = pd.read_csv(Ldir['run'] + 'rname_list.txt', header=None,
        names=['River Name'])
    rnames = rdf['River Name'].values
    rnames = rnames.tolist()
    rnames[rnames.index('duwamish')] = 'green'
    rnames[rnames.index('hammahamma')] = 'hamma'
    rnames.remove('skagit_south')
else:
    # override for a shorter list
    rnames = ['dosewallips']

from datetime import datetime
dt0 = datetime(1980,1,1)
dt1 = datetime(2015,12,31)
days = (dt0, dt1)

for rn in rnames:
   
    riv = river_class.River(Ldir)
    riv.name_it(rn)   
    riv.get_ecy_info(riv.name)   
    riv.get_nws_info(riv.name)
    if rn == 'fraser':
        riv.get_ec_data(days)
    else:
        riv.get_usgs_data(days)
    riv.print_info()
    
    rqt = riv.qt
    rqt.to_pickle(outdir + 'clim_' + rn + '.p')
    

    




