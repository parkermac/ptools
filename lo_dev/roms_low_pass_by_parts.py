"""
Makes a low-passed version of a series of ROMS history files.
"""

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
import time
import numpy as np

# create the filter
filt0 = zfun.godin_shape()
nf = len(filt0)

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
            
ds = nc.MFDataset(flist)
vname = 'salt'

# TIMING #
tt0 = time.time()

# define a function that does part of the filtering
def get_part(i0, i1, ds, vname, filt0):     
    v = ds.variables[vname][i0:i1,:,:,:]
    print('\n' + 30*'*')    
    print('i0 = ' + str(i0))
    print('i1 = ' + str(i1))
    print('v shape = ' + str(v.shape))
    sys.stdout.flush()
    filt1 = filt0[i0:i1]
    filt = filt1.reshape((len(filt1),1,1,1))
    vf = (filt*v).sum(axis=0)
    return vf

# set a number of time blocks to use
nblocks = 1 # works perfectly for 1 to 71

nb = int(nf/nblocks)

# initial values
i0 = 0
i1 = nb
vf = get_part(i0, i1, ds, vname, filt0)
# loop over the rest
if nblocks > 1:
    while i1 < nf:        
        i0 += nb
        i1 += nb    
        if i1 >= nf+1:
            i1 = nf+1       
        vf += get_part(i0, i1, ds, vname, filt0)
ds.close()

# TIMING #
print('\n' + 30*'=')
print('Took ' + str(round(time.time() - tt0)) + ' seconds')
print('Mean vf = ' + str(vf.mean()))
sys.stdout.flush()
