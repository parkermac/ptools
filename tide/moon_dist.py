"""
Code to test ephem moon distance.

"""

import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import pytz
import ephem
import pandas as pd

moon = ephem.Moon()

tz_utc = pytz.timezone('UTC')
dt = datetime(2016,1,1,tzinfo=tz_utc)

dt_list = []
md_list = []

while dt <= datetime(2016,12,31,tzinfo=tz_utc):


    moon.compute(ephem.Date(dt))
    md = moon.earth_distance
    dt_list.append(dt)
    md_list.append(md)
    
    dt = dt + timedelta(days=1)

mds = pd.Series(index=dt_list, data=md_list)
# plotting
plt.close('all')

mds.plot()

plt.show()

