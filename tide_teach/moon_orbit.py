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
dt1 = datetime(year+19,12,31,tzinfo=pytz.timezone('UTC'))

def add_tf(df):
    # add the (normalized) tractive force
    r = 1000 * df['Distance (km)'].values # distance in m
    tf = tfun.get_tractive_scale(r)
    tfn = tf/tf.mean() # normalize tractive force by its mean
    df['Tractive Force'] = tfn
    return df

moon_orbit_df = efun.get_moon_orbit(dt0, dt1)
moon_robit_df = add_tf(moon_orbit_df)

fm_df, nm_df = efun.get_full_new(dt0, dt1)

# plotting
plt.close('all')

moon_orbit_df.loc[:,['Declination (deg)', 'Tractive Force']].plot(subplots=True, grid=True)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# fm_df.plot(y='Distance (km)', label='Full Moon Distance (km)', style='-or', ax=ax)
# nm_df.plot(y='Distance (km)', label='New Moon Distance (km)', style='-ob', ax=ax)

plt.show()

