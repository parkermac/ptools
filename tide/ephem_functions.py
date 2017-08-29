#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 13:29:04 2017

@author: pm7

Functions that make use of eht ephem module.

"""
from datetime import datetime
import pytz
import ephem

def make_info(city='Seattle', zone='US/Pacific'):
    # get timezone and observer info
    # use pytz.common_timezones to see a list of possibilities
    tz_local = pytz.timezone(zone)
    tz_utc = pytz.timezone('UTC')
    if city == 'Westport':
        # make obs by hand as a modification of Seattle
        obs = ephem.city('Seattle')
        obs.lon = '-124:06:36'
        obs.lat = '46:53:27'
    else:
        # or look if up in the ephem database (e.g. for Seattle)
        obs = ephem.city(city)
    return tz_utc, tz_local, obs

def get_sunmoon(dt_local, tz_utc, obs):
    # create sun and moon objects referenced to an observer
    # and a local datetime
    dt_utc = dt_local.astimezone(tz_utc)
    eph_date = ephem.Date(dt_utc)
    obs.date = eph_date
    sun = ephem.Sun(obs)
    moon = ephem.Moon(obs)
    return sun, moon
    
def epht_to_local(epht, tz_local):
    # convert an ephem time to local datetime
    dt_naive = epht.datetime()
    dt_utc = pytz.utc.localize(dt_naive)
    dt_local = dt_utc.astimezone(tz_local)
    return dt_local
    
def get_times(dt_local, tz_utc, tz_local, obs):
    sun, moon = get_sunmoon(dt_local, tz_utc, obs)
    S = dict()
    M = dict()
    S['rise'] = epht_to_local(obs.next_rising(sun), tz_local)
    S['transit'] = epht_to_local(obs.next_transit(sun), tz_local)
    S['set'] = epht_to_local(obs.next_setting(sun), tz_local)
    M['rise'] = epht_to_local(obs.next_rising(moon), tz_local)
    M['transit'] = epht_to_local(obs.next_transit(moon), tz_local)
    M['set'] = epht_to_local(obs.next_setting(moon), tz_local)
    return S, M
    
def print_times(dt_local, S, M, tz_local):
    # print results
    fmt = '%14s %25s %s'
    tlist = ['rise', 'transit', 'set']
    print(fmt % ('Date', dt_local.ctime(), tz_local.zone))
    for tt in tlist:
        print(fmt % (('* sun '+tt), S[tt].ctime(), tz_local.zone))
    for tt in tlist:
        print(fmt % (('o moon '+tt), M[tt].ctime(), tz_local.zone))
          
if __name__ == '__main__':
    # example use of the functions
    tz_utc, tz_local, obs = make_info()
    dt_local = datetime(2016, 9, 19, tzinfo=tz_local)
    S, M = get_times(dt_local, tz_utc, tz_local, obs)
    print_times(dt_local, S, M, tz_local)
