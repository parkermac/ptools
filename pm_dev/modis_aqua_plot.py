"""
Plot MODIS Aqua Data

"""
import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()
import zfun

import netCDF4 as nc

indir = Ldir['parent'] + 'ptools_data/modis_aqua/'
#fn = 'A20172332017240.L3m_8D_CHL_chl_ocx_9km.nc'
fn = 'A20172332017240.L3b_8D_CHL.nc' # can't get a proper listing...

ds = nc.Dataset(indir + fn)

zfun.ncd(ds)