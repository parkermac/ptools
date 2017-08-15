#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 14:36:08 2017

@author: PM5

Test of the new astral and ephem modules.
"""
from datetime import datetime
import pytz
import ephem

# get timezone info
# use pytz.common_timezones to see a list of possibilities
tz_pac = pytz.timezone('US/Pacific')
tz_utc = pytz.timezone('UTC')

# set an initial datetime
dt_utc0 = datetime(2016,1,1,tzinfo=tz_utc)

# and convert it to an ephem Date object
eph_date0 = ephem.Date(dt_utc0)
# Run str() on this to see the UTC date it represents.

# create the ephem Observer object
city = 'Seattle'
obs = ephem.city(city)
# and use these to get things like rise and set time
sun = ephem.Sun(obs)
moon = ephem.Moon(obs)

def get_next_full_moon(eph_date):
    eph_date = ephem.next_full_moon(eph_date)
    return eph_date
    
def print_info(city, eph_date):
    dt_naive = eph_date.datetime()
    # although eph dates are always UTC, when you first convert them to
    # datetime they are timezone naive
    dt_utc = pytz.utc.localize(dt_naive) # a key step
    dt_pac = dt_utc.astimezone(tz_pac)    
    print(city + ' next full moon: ' + dt_utc.ctime() + ' UTC')
    print(city + ' next full moon: ' + dt_pac.ctime() + ' Pacific')    
    return dt_utc

# output
print('\nStarting time: ' + dt_utc0.ctime() + ' UTC')

# initialize 
eph_date = get_next_full_moon(eph_date0)
ii = 0
# loop over a year
while eph_date.datetime().year < 2017:
    # print info
    print('\n** ii = ' + str(ii))
    dt_utc = print_info(city, eph_date)
    
    # this finds the sunrise time on that day
    # or the day after
    dt_floor = datetime(dt_utc.year, dt_utc.month, dt_utc.day)
    eph_date_floor = ephem.Date(dt_floor)
    obs.date = eph_date_floor
    sun = ephem.Sun(obs)
    rt = sun.rise_time
    drt_naive = rt.datetime()
    drt_utc = pytz.utc.localize(drt_naive)
    drt_pac = drt_utc.astimezone(tz_pac)    
    print('Sunrise at: ' + drt_pac.ctime() + ' Pacific')
    
    # increment eph_date
    eph_date = get_next_full_moon(eph_date)
    ii += 1
    



