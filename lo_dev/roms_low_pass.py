"""
Makes a low-passed version of a series of ROMS history files.

RESULT: it filters 71 files (phys variables, 8 are 4D) in ~1 minute.
    => 6 hours for a year.
"""

# TIMING #
import time
tt0 = time.time()

gridname = 'cascadia1'
tag = 'base'
ex_name = 'lo1'

# setup
import os; import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(gridname, tag)
Ldir['gtagex'] = Ldir['gtag'] + '_' + ex_name
import netCDF4 as nc
import zfun; reload(zfun)

# create the list of history files
date_string = '2015.09.19'
indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + date_string + '/')
if False:
    # You can call MFDataset with a wild card 
    flist = indir + 'ocean_his_*.nc'
else:
    # or you can call it with a list of files
    # full length for godin would be (2,73) = 71 elements
    flist = []
    for ii in range(2,73): 
        hnum = ('0000' + str(ii))[-4:]
        flist.append(indir + 'ocean_his_' + hnum + '.nc')
       
# create the filter
nf = len(flist)
if nf == 71:
    print('Godin filter')
    filt0 = zfun.godin_shape()
else:
    print('Hanning filter for list length = ' + str(nf))
    filt0 = zfun.hanning_shape(nf)
#    
# create the output file
fnout = '../../ptools_output/low_pass_test.nc'
import shutil
shutil.copyfile(flist[0],fnout)
#
# create the Datasets                
ds = nc.MFDataset(flist)
dsout = nc.Dataset(fnout)
#
# loop over all variables that have time axes
for vn in ds.variables:
    if 'ocean_time' in ds.variables[vn].dimensions:
        print(vn + ' ' + str(ds.variables[vn].shape))
        ndim = len(ds.variables[vn].shape)
        filt_shape = (nf,)
        for ii in range(ndim-1):
            filt_shape = filt_shape + (1,)
        v = ds.variables[vn][:]
        filt = filt0.reshape(filt_shape)
        vf = (filt*v).sum(axis=0)
        dsout.variables[vn] = vf
ds.close()
dsout.close()

# TIMING #
print('\n' + 30*'=')
print('Took ' + str(round(time.time() - tt0)) + ' seconds')
sys.stdout.flush()
