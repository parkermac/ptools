"""
Code to experiment with getting HYCOM data from a forecast.

Parker MacCready
"""

# setup
import os; import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart(alp)

alpo = os.path.abspath('../../LiveOcean/setup/make_forcing/ocn')
if alpo not in sys.path:
    sys.path.append(alpo)
import Ofun; reload(Ofun)

print '** START getting catalog'
# create a list of url's of the preprocessed HYCOM files for this forecast
fn_list = Ofun.get_hycom_file_list()
print '** DONE getting catalog'

# get a selection of the raw list (e.g. one per day)
var_list = ['ssh','ts3z','uv3z']
varf_dict = dict()
for var_name in var_list:
    varf_dict[var_name] = Ofun.make_shortened_list(fn_list, var_name)
# check that all the lists are the same length
list_len = []
for var_name in var_list:
    list_len.append(len(varf_dict[var_name]))
if list_len.count(list_len[0]) != len(list_len):
    print 'WARNING: Lists in varf_dict are different lengths!'
else:
    NT = list_len[0] # the number of times, to use below
# really we should be checking that the times in each list item are identical
    
# for each variable type download a selected region
#
# this took about 30 sec to get all variables at one time, so I expect it
# will take 5-6 minutes to get all ~11 times, but performance has been spotty
# today - mainly getting the catalog xml...

# define the output location
nc_dir = '../../tools_output/pydev_out/hycom_test/'

# wow this is very cool (but OBSOLETE)
nt = 0
while nt < NT:
    out_dict_allvar = dict()
    for var_name in var_list:
        fn = varf_dict[var_name][:][nt]
        out_dict_allvar[var_name] = Ofun.get_extraction(fn, var_name)
    if nt == 0:
        Ofun.initialize_netcdf_files(out_dict_allvar, nc_dir)
    else:
        Ofun.append_netcdf_files(out_dict_allvar, nt, nc_dir)
    nt += 1

import netCDF4 as nc
for vn in ['ssh', 's3d', 't3d', 'u3d', 'v3d']:
    print '*************** ' + vn + '*******************************'
    ds = nc.Dataset(nc_dir + vn + '.nc')
    for vn in ds.variables:
        print ds.variables[vn]
    ds.close()
       



