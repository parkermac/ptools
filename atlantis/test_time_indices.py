# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:05:40 2016

@author: PM5

Get info on timing and indices of low-passed files.

"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun
import zrfun
import collections

from datetime import datetime
import matplotlib.dates as mdates

whichyear = 2005

if Ldir['env'] == 'pm_mac': # mac version
    in_dir = Ldir['parent'] + 'roms/output/salish_2006_4/'
elif Ldir['env'] == 'fjord': # fjord version
    if whichyear == 2006:
        in_dir = '/pmr3/pmraid1/daves/runs/salish_2006_4/OUT/'
        out_dir0 = '/boildat1/parker/roms/output/salish_2006_4_lp/'
    elif whichyear == 2005:
        in_dir = '/pmr3/pmraid2/daves/salish_2005_1/OUT/'
        out_dir0 = '/boildat1/parker/roms/output/salish_2005_1_lp/'
        
for ii0 in range(26, 8666 + 1, 24):
#for ii0 in range(26, 300 + 1, 24):
    start_time = datetime.now()

    NF = 71
    flist = []
    for ii in range(NF):
        hh = ('000' + str(ii0 + ii))[-4:]
        flist.append(in_dir + 'ocean_his_' + hh + '.nc')

    tm_list = []
    for fn in flist:
        T = zrfun.get_basic_info(fn, only_T=True)
        tm_list.append(T['tm'])
        
    # find the mean time
    mdv = mdates.date2num(tm_list)
    TM = mdates.num2date(mdv.mean())

    # make output name (full path) using LiveOcean naming convention
    print('%s %d' % (TM.strftime('%Y.%m.%d'), ii0) )