"""
Code to find local times of sun and moon
rise, transit, and set on a given day.

"""

from datetime import datetime, timedelta
import pytz
import ephem

def eph_to_local(epht, tz_local):
    dt_naive = epht.datetime()
    dt_utc = pytz.utc.localize(dt_naive)
    dt_local = dt_utc.astimezone(tz_local)
    return dt_local

def get_sunmoon(dt_local, tz_utc, obs):
    dt_utc = dt_local.astimezone(tz_utc)
    eph_date = ephem.Date(dt_utc)
    obs.date = eph_date
    sun = ephem.Sun(obs)
    moon = ephem.Moon(obs)
    return sun, moon
    
def get_dt0(dt, tz):
    dt0 = datetime(dt.year, dt.month, dt.day, tzinfo=tz)
    return dt0
    
def get_times(dt_local, tz_utc, tz_local, obs):
    sun, moon = get_sunmoon(dt_local, tz_utc, obs)
    S = dict()
    M = dict()
    S['rise'] = eph_to_local(obs.next_rising(sun), tz_local)
    S['transit'] = eph_to_local(obs.next_transit(sun), tz_local)
    S['set'] = eph_to_local(obs.next_setting(sun), tz_local)
    M['rise'] = eph_to_local(obs.next_rising(moon), tz_local)
    M['transit'] = eph_to_local(obs.next_transit(moon), tz_local)
    M['set'] = eph_to_local(obs.next_setting(moon), tz_local)
    return S, M
    
def print_times(dt_local, S, M, tz_local):
    # print results
    fmt = '%14s %25s %s'
    tlist = ['rise', 'transit', 'set']
    print(fmt % ('Date', dt_local.ctime(), tz_local.zone))
    for tt in tlist:
        print(fmt % (('-sun '+tt), S[tt].ctime(), tz_local.zone))
    for tt in tlist:
        print(fmt % (('-moon '+tt), M[tt].ctime(), tz_local.zone))

def make_info(city='Seattle', zone='US/Pacific'):
    # get timezone info
    # use pytz.common_timezones to see a list of possibilities
    tz_local = pytz.timezone(zone)
    tz_utc = pytz.timezone('UTC')
    # create the ephem Observer object
    obs = ephem.cities.lookup(city) # check for mistakes!
    # or use this version for cities in the ephem database
    # obs = ephem.city(city)
    return tz_utc, tz_local, obs
        
def make_dt_local(y, m, d, tz_local):
    # get local datetime
    dt_local = datetime(y, m, d, tzinfo=tz_local)
    return dt_local
    
if __name__ == '__main__':
    tz_utc, tz_local, obs = make_info()
    dt_local = make_dt_local(2016, 6, 20, tz_local)
    S, M = get_times(dt_local, tz_utc, tz_local, obs)
    print_times(dt_local, S, M, tz_local)

