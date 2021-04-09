"""
Create a tidally averaged surface Bernoulli function
from a LiveOcean layer extraction.
"""

# setup
import sys
from pathlib import Path
import netCDF4 as nc
import numpy as np

pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
Ldir = Lfun.Lstart()
import plotting_functions as pfun

in_fn = Ldir['parent'] / 'LiveOcean_output' / 'layer' / 'cas6_v3_lo8b_2019.06.01_2019.08.31' / 'surface_hourly.nc'
# a total of 92 days, so we can form 90 tidal averages I believe

ds = nc.Dataset(in_fn)

for vn in ds.variables:
    print('%s %s' % (vn, str(ds[vn].shape)))

