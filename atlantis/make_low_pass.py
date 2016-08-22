# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:05:40 2016

@author: PM5

Code to make low-passed files for the Atlantis project.

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

from datetime import datetime, timedelta
start_time = datetime.now()

in_dir = Ldir['parent'] + 'roms/output/salish_2006_4/'
out_dir0 = Ldir['parent'] + 'roms/output/salish_2006_4_lp/'
Lfun.make_dir(out_dir0)

# make input list (full paths)
flist = []
# create the list of history files
#
# for salish_2006_4 we have 0025-8715, which have times:
# 0025: 'tm': datetime.datetime(2006, 1, 3, 0, 0)
# 4993: 'tm': datetime.datetime(2006, 7, 29, 0, 0)
# 8715: 'tm': datetime.datetime(2006, 12, 31, 2, 0)
#
# starting at 4994 (1:00 AM, 2006.07.29) with NF = 71 then returns a file with
# 'tm': datetime.datetime(2006, 7, 30, 12, 0) = GOOD
#
# for the full year we should start at ii0 = 26 and go to 8666, like this:
# a = range(26, 8666 + 1, 24)

for ii0 in [4994]:

    NF = 71
    for ii in range(NF):
        hh = ('000' + str(ii0 + ii))[-4:]
        flist.append(in_dir + 'ocean_his_' + hh + '.nc')

    tm_list = []
    for fn in flist:
        [T] = zrfun.get_basic_info(fn, getG=False, getS=False)
        tm_list.append(T['tm'])
    # find the mean time
    import matplotlib.dates as mdates
    mdv = mdates.date2num(tm_list)
    TM = mdates.num2date(mdv.mean())

    # make output name (full path) using LiveOcean naming convention
    out_dir = out_dir0 + 'f' + TM.strftime('%Y.%m.%d') + '/'
    Lfun.make_dir(out_dir)
    out_fn = (out_dir + '/low_passed.nc')

    # create the filter
    nf = len(flist)
    if nf == 71:
        print(' - Using Godin filter')
        filt0 = zfun.godin_shape()
    else:
        print(' - Using Hanning filter for list length = ' + str(nf))
        filt0 = zfun.hanning_shape(nf)

    # RUN THE FUNCTION
    zrfun.roms_low_pass(flist, out_fn, filt0)

    # save result info
    result_dict = collections.OrderedDict()
    time_format = '%Y.%m.%d %H:%M:%S'
    result_dict['start_time'] = start_time.strftime(time_format)
    end_time = datetime.now()
    result_dict['end_time'] = end_time.strftime(time_format)
    dt_sec = (end_time - start_time).seconds
    result_dict['total_seconds'] = str(dt_sec)
    if os.path.isfile(out_fn):
        result_dict['result'] = 'success'
    else:
        result_dict['result'] = 'fail'
    info_dir = out_dir + 'Info/'
    Lfun.make_dir(info_dir)
    Lfun.dict_to_csv(result_dict, info_dir + 'process_status.csv')
