"""
Code to download WCOFS model output, for the COMT 3 project.

This is for the newer version, accessed March 9, 2020...
"""
import requests
from time import time

tt0 = time()
fn =  '/thredds/fileServer/NOAA/WCOFS/MODELS/201901/nos.wcofs.avg.nowcast.20190101.t03z.nc'
url = 'http://opendap.co-ops.nos.noaa.gov' + fn

local_filename = url.split('/')[-1]
print(local_filename)
# NOTE the stream=True parameter below

with requests.get(url, stream=True) as r:
    r.raise_for_status()
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192): 
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
                
print('took %0.1f seconds' % (time()-tt0))

