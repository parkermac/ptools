"""
Code to plot tides from NOAA csv files.

It makes a MOVIE where each frame is the tide height over a day,
so you can see how low tide marches forward avery day.

"""
import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

from importlib import reload
import ephem_functions as efun
reload(efun)

# User choices
season = 'summer' # summer, winter
city = 'Seattle' # Seattle, Westport
testing = False

home = os.environ.get('HOME')
dir00 = home + '/Documents/'

# read in tide data (time is UTC)
indir = dir00 + 'ptools_data/tide/'
if city == 'Seattle':
    fn = 'CO-OPS__9447130__hr.csv'
    city = 'Seattle'
    zone='US/Pacific'
elif city == 'Westport':
    fn = 'CO-OPS__9441102__hr.csv'
    city = 'Westport'
    zone='US/Pacific'
    
# Date Time, Water Level, Sigma, I, L
# 2016-01-01 00:00,1.302,0.000,0,0
# time is UTC
df = pd.read_csv(indir+fn, index_col='Date Time', parse_dates = True)
for k in df.keys():
    df = df.rename(columns={k: k.strip()})
df = df.drop(['Sigma', 'I', 'L'], axis=1)
df = df.rename(columns={'Water Level': 'SSH_obs'})

# PLOTTING LOOP

# initializing
plt.close('all')
if season == 'summer':
    d0 = datetime(2016,6,10,7,0,0) # Summer start (PDT)
    tag = season
elif season == 'winter':
    d0 = datetime(2016,11,30,8,0,0) # Winter start (PST)
    tag = season
# starting on hour 8 of UTC is midnight PST
d1 = d0 + timedelta(days=1)
d2 = d1 + timedelta(days=1)
# we get data from two days in order to do the blending

# set the number of days to plot
if testing == False:
    # 29.57 days is two S-N cycles
    d_last = d0 + timedelta(days=29)
elif testing == True:
    # testing
    d_last = d0 + timedelta(days=5)
    
# prepare a directory for results
outdir0 = dir00 + 'ptools_output/tide/'
Lfun.make_dir(outdir0, clean=False)
outdir = outdir0 + 'day_movie_' + city + '_' + tag + '/'
Lfun.make_dir(outdir, clean=True)

# RC SETUP (plotting defaults)
def set_rc(fs, lw, mks):
    fs_big = fs
    fs_small = fs-4
    lw_big = lw
    lw_small = lw - 2
    plt.rc('xtick', labelsize=fs_small)
    plt.rc('ytick', labelsize=fs_small)
    plt.rc('xtick.major', size=10, pad=5, width=lw_small)
    plt.rc('ytick.major', size=10, pad=5, width=lw_small)
    plt.rc('axes', lw=lw_small)
    plt.rc('lines', lw=lw_big, markersize=mks)
    plt.rc('font', size=fs_big)
    plt.rc('grid', color='g', ls='-', lw=lw_small, alpha=.3)
fs = 20
lw = 5
mks = 50
set_rc(fs, lw, mks)

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1,1,1)
ax.grid()
ax.set_xlabel('Local Time of Day')
ax.set_ylabel('Tide Height (m)')
ax.set_xlim(0, 24)
ax.set_xticks([6, 12, 18])
ax.set_xticklabels(['6AM', 'Noon', '6PM'])
y0 = -1
y1 = 5
ax.set_ylim(y0, y1)

tz_utc, tz_local, obs = efun.make_info(city=city, zone=zone)

i_plot = 0 # initialize plot name number
while d1 <= d_last:
    # get sun and moon info
    dt0_local = datetime(d0.year, d0.month, d0.day, tzinfo=tz_local)
    S0, M0 = efun.get_times(dt0_local, tz_utc, tz_local, obs)
    mth0 = M0['transit'].hour + M0['transit'].minute/60
    srh0 = S0['rise'].hour + S0['rise'].minute/60
    ssh0 = S0['set'].hour + S0['set'].minute/60
    dt1_local = datetime(d1.year, d1.month, d1.day, tzinfo=tz_local)
    S1, M1 = efun.get_times(dt1_local, tz_utc, tz_local, obs)
    mth1 = M1['transit'].hour + M1['transit'].minute/60
    srh1 = S1['rise'].hour + S1['rise'].minute/60
    ssh1 = S1['set'].hour + S1['set'].minute/60
    # avoid blending when we jump from hour 24 to hour 0
    if mth1 < mth0:
        mth0 = np.nan
        mth1 = np.nan
    
    # get one-day chunks of the tide height time series
    a = df.loc[d0:d1].values
    b = df.loc[d1:d2].values
    aa = np.array(a).flatten()
    bb = np.array(b).flatten()
    tt = np.arange(0,len(aa))
    
    ndiv = 8 # number of divisions of the day (for a smoother movie)
    i_fr = 0 # index into the divisions of the day
    while i_fr < ndiv:
        # name the output file
        nouts = ('0000' + str(i_plot))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir + outname
        if (i_plot/10).is_integer():
            print(outname)
        # calculate fraction of the way through the day for this fames
        fr = i_fr/ndiv
        # tide signal at this time - blended to make movie smooth
        cc = (1-fr)*aa + fr*bb
        # hours of sun rise and set - blended
        srh = (1-fr)*srh0 + fr*srh1
        ssh = (1-fr)*ssh0 + fr*ssh1
        # hour of moon efun - blended
        mth = (1-fr)*mth0 + fr*mth1
        # plot various things
        lh_night0 = ax.fill([0, srh, srh, 0], [y0, y0, y1, y1], 'k', alpha=.3)
        lh_night1 = ax.fill([ssh, 24, 24, ssh], [y0, y0, y1, y1], 'k', alpha=.3)
        lh_day = ax.fill([srh, ssh, ssh, srh], [y0, y0, y1, y1], 'y', alpha=.3)
        lh_moon = ax.plot(mth, y1-1, 'o', color='w', markeredgecolor='k')
        lh_tide = ax.plot(tt, cc, '-b')
        lh_minus = ax.fill_between(tt, cc, where=(cc<=0), interpolate=True, color='red')
        lh_date = ax.text(2, y1-.5, city + d0.strftime(': %B %d'), weight='bold')
        plt.tight_layout()
        plt.draw()
        plt.savefig(outfile)
        # remove the line and other things
        lh_night0.pop(0).remove()
        lh_night1.pop(0).remove()
        lh_day.pop(0).remove()
        lh_moon.pop(0).remove()
        lh_tide.pop(0).remove()
        lh_minus.remove()
        lh_date.remove()
        # increment counters
        i_fr += 1
        i_plot += 1
    # increment day
    d0 = d1
    d1 = d0 + timedelta(days=1)
    d2 = d1 + timedelta(days=1)

# RC CLEANUP
plt.rcdefaults()

# and make a movie (does not work in Spyder)
import subprocess
cmd = ['ffmpeg',
    '-r', '24', # framerate fps
    '-i', outdir+'plot_%04d.png',
    '-vcodec', 'libx264',
    '-pix_fmt', 'yuv420p',
    '-crf', '25',
    outdir + 'movie.mp4']
proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if False:
    print('\n-main: screen output from subprocess-')
    print(proc.stdout.decode())
    print('\n-main: errors from subprocess-')
    # for some reason the ffmpeg output ends up in stderr
    print(proc.stderr.decode())

