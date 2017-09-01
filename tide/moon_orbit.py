"""
Code to explore moon orbital properties.

"""

import matplotlib.pyplot as plt
from datetime import datetime, timedelta, timezone
import pandas as pd
import pytz

from importlib import reload
import ephem_functions as efun
reload(efun)
import tractive_functions as tfun
reload(tfun)

year = 2016
dt0 = datetime(year,1,1,tzinfo=pytz.timezone('UTC'))
dt1 = datetime(year,12,31,tzinfo=pytz.timezone('UTC'))

moon_orbit_df = efun.get_moon_orbit(dt0, dt1)

r = 1000 * moon_orbit_df['Distance (km)'].values # distance in m
tf = tfun.get_tractive_scale(r)
tfp = 100 * (tf/tf.mean() - 1)

moon_orbit_df['Tractive Force (% from mean)'] = tfp

fm_df, nm_df = efun.get_full_new(dt0, dt1)
a = pd.concat([fm_df, nm_df])
a = a.sort_index()

# plotting
plt.close('all')

moon_orbit_df.loc[:,['Declination (deg)', 'Tractive Force (% from mean)']].plot(grid=True)

if True:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    a['Distance (km)'][a['fullnew']==1].plot(style='-or', ax=ax)
    a['Distance (km)'][a['fullnew']==0].plot(style='-ob', ax=ax)

plt.show()

