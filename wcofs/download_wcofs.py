"""
Code to download WCOFS model output, for the COMT 3 project.

This is for the newer version, accessed March 9, 2020...

RESULT: it works but performance is terrible.  Almost 10 minutes on both boiler and my mac
for a file that is just 361 MB.  Thypcially this would be 36 sec at home and 4 sec at work.
"""

import os; import sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun

import requests
from datetime import datetime, timedelta
from time import time

def boolean_string(s):
    if s not in ['False', 'True']:
        raise ValueError('Not a valid boolean string')
    return s == 'True' # note use of ==

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# standard arguments
parser.add_argument('-0', '--date_string0', type=str, default='2019.01.01')
parser.add_argument('-1', '--date_string1', type=str, default='2019.01.03')
parser.add_argument('-test', '--testing', type=boolean_string, default=True)
args = parser.parse_args()

# Get Ldir
Ldir = Lfun.Lstart()

# create the output directory
outdir0 = Ldir['roms'] + 'output/'
Lfun.make_dir(outdir0)
outdir = outdir0 + 'wcofs_avg_now/'
Lfun.make_dir(outdir)

# get time limits
ds0 = args.date_string0; ds1 = args.date_string1
dt0 = datetime.strptime(ds0, '%Y.%m.%d')
dt1 = datetime.strptime(ds1, '%Y.%m.%d')

# create list of days (datetimes)
dt_list = []
dt = dt0
while dt <= dt1:
    dt_list.append(dt)
    dt += timedelta(days=1)

for dt in dt_list:
    this_mo = datetime.strftime(dt, '%Y%m')
    this_day = datetime.strftime(dt, '%Y%m%d')
    print('\nWorking on: %s' % (this_day))
    
    tt0 = time()
    fn =  ('/thredds/fileServer/NOAA/WCOFS/MODELS/' + this_mo +
        '/nos.wcofs.avg.nowcast.' + this_day + '.t03z.nc')
    url = 'http://opendap.co-ops.nos.noaa.gov' + fn
    outfile = outdir + url.split('/')[-1]
    print(' ' + url)
    print(' ' + outfile)

    if args.testing == False:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(outfile, 'xb') as f:
                for chunk in r.iter_content(chunk_size=8192): 
                    if chunk: # filter out keep-alive new chunks
                        f.write(chunk)
            
    print('  -- took %0.1f seconds' % (time()-tt0))
    sys.stdout.flush()

