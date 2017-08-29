"""
Code to test ephem moon distance.

"""

import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pytz
import ephem
import pandas as pd
import numpy as np

moon = ephem.Moon()

tz_utc = pytz.timezone('UTC')
year = 2017

dt0 = datetime(year,1,1,tzinfo=tz_utc)
dt1 = datetime(year,12,31,tzinfo=tz_utc)

dt_list = []
dist_list = []
phase_list = []
dec_list = []

dt = dt0
while dt <= dt1:
    
    moon.compute(ephem.Date(dt))
    
    dt_list.append(dt)
    
    dist_list.append(moon.earth_distance)
    # distance in AU = astronomical units = 149597871 km
    
    phase_list.append(moon.moon_phase)
    # phase as fraction of the surface illuminated
    
    dec_list.append(moon.dec)
    # declination in radians
    
    dt = dt + timedelta(days=1)

ddict = {'Dist':dist_list, 'Ph':phase_list, 'Dec':dec_list}
mdf = pd.DataFrame(index=dt_list, data=ddict)

# also get new and full moons
fm_list = []
nm_list = []

dt = dt0
while dt <= dt1:
    eph_date = ephem.Date(dt)
    eph_date = ephem.next_full_moon(eph_date)
    dt_naive = eph_date.datetime()
    dt = pytz.utc.localize(dt_naive)
    if dt <= dt1:
        fm_list.append(dt)
    
dt = dt0
while dt <= dt1:
    eph_date = ephem.Date(dt)
    eph_date = ephem.next_new_moon(eph_date)
    dt_naive = eph_date.datetime()
    dt = pytz.utc.localize(dt_naive)
    if dt <= dt1:
        nm_list.append(dt)

fm = pd.DataFrame(index=fm_list, data={'Full':np.ones(len(fm_list))})
nm = pd.DataFrame(index=nm_list, data={'New':np.zeros(len(nm_list))})

a = pd.concat([mdf, fm, nm])
a = a.sort_index()

# plotting
plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(111)
a['Full'].plot(style='*r', ax=ax)
a['New'].plot(style='ob', ax=ax)
a['Ph'].plot(style='-k', ax=ax)

plt.show()

