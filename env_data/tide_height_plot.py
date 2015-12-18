"""
This plots tide height data.
"""

dir0 = '/Users/PM5/Documents/'

# setup
import os; import sys
alp = os.path.abspath(dir0 + 'LiveOcean/alpha/')
if alp not in sys.path:
    sys.path.append(alp)
import zfun; reload(zfun)
import matfun; reload(matfun)
    
indir = dir0 + 'tools_data/obs_data/tide_gauge/NOAA/NOAA_2013/'

dl = os.listdir(indir)
dl2 = []
for item in dl:
    if '.mat' in item:
        #print(item)
        dl2.append(item)
        
for a in dl2:
    aa = matfun.loadmat(indir + a)['tgauge']
    print(a + ': ' + aa['station_name'])
        

