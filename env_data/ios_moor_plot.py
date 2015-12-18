"""
Code to look at the contents of the IOS mooring records.
"""
# setup
import os; import sys
alp = os.path.abspath('../../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun; reload(zfun)
import netCDF4 as nc
from datetime import datetime, timedelta

dir0 = '/Users/PM5/Documents/tools_data/obs_data/mooring/IOS/'

def matlabdn2datetime(dn):
    dt = (datetime.fromordinal(int(dn)) +
        timedelta(days=dn%1) - timedelta(days = 366))
    return dt
        
def datetime2matlabdn(dt):
   mdn = dt + timedelta(days = 366)
   frac_seconds = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   frac_microseconds = dt.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
   return mdn.toordinal() + frac_seconds + frac_microseconds

for mnum in range(1,9):
    
    fn = dir0 + 'IOS_adcp' + str(mnum) + '.nc'
    
    ds = nc.Dataset(fn)
    
    #zfun.ncd(ds)
    
    vn_list = ['pressure','temperature','salinity',
        'current_velocity_along_shore','current_velocity_across_shore',
        'time','latitude','longitude']
    
    aa = dict()   
    for vn in vn_list:
        aa[vn] = ds.variables[vn][:]
    ds.close()
    
    print 50*'*'
    print fn    
    print 'longitude = ' + str(aa['longitude'][0,0])
    print 'latitude = ' + str(aa['latitude'][0,0])
    for ii in range(len(aa['pressure'])):
        print 'pressure = ' + str(aa['pressure'][ii])
    
    
    # note on time units
    # "serial decimal starting Jan 1 0000"

    dt_list = []
    for ii in range(len(aa['time'])):
        dt_list.append(matlabdn2datetime(aa['time'][ii]))
    
    dn0 = aa['time'][0]
    dn1 = aa['time'][-1]    
    dt0 = dt_list[0]
    dt1 = dt_list[-1]
    print 'start time = ' + dt0.strftime('%Y.%m.%d %H:%M:%S')
    print 'end time = ' + dt1.strftime('%Y.%m.%d  %H:%M:%S')





