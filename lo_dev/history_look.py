"""
Code for exploring ROMS history files.
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
#import netCDF4 as nc
import zfun; reload(zfun)

# create the list of history files
date_string = '2015.09.19'
indir = (Ldir['roms'] + 'output/' + Ldir['gtagex'] + '/f' + date_string + '/')
fn = indir + 'ocean_his_0002.nc'

G, S, T = zfun.get_basic_info(fn)

nr, nc = G['lon_rho'].shape

print('\n' + fn)
print('nr = %d, nc = %d' % (nr, nc))