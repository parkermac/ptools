"""
Code to test the timespan covered by EC river data sources, and attempts at error handling.
"""

import requests
import bs4
import pandas as pd
import sys

from datetime import datetime

# functions lightly modified from LiveOcean/forcing/riv2/river_class.py

def get_ec_data(days, timeout=10):
    # gets Environment Canada data, using code cribbed from:
    #https://bitbucket.org/douglatornell/ecget/src/
    # NOTE: this will get data up through today, but can only go back
    # 18 months into the past.
    # To get longer records use get_ec_data_historical.
    import requests
    import bs4
    try:
        PARAM_IDS = {'discharge': 47,}

        params = {
            'mode': 'Table',
            'type': 'realTime',
            'prm1': PARAM_IDS['discharge'],
            'prm2': -1,
            'stn': '08MF005',
            'startDate': days[0].strftime('%Y-%m-%d'),
            'endDate': days[1].strftime('%Y-%m-%d'),
        }
        DATA_URL = 'http://wateroffice.ec.gc.ca/report/real_time_e.html'
        DISCLAIMER_COOKIE = {'disclaimer': 'agree'}
        response = requests.get(DATA_URL, params=params,
                                cookies=DISCLAIMER_COOKIE)
        soup = bs4.BeautifulSoup(response.content, 'lxml')
        table = soup.find('table')
        table_body = table.find('tbody')
        rows = table_body.find_all('tr')
        data = []
        for row in rows:
            cols = row.find_all('td')
            cols = [ele.text.strip() for ele in cols]
            data.append([ele for ele in cols if ele]) # Get rid of empty values
        d_dict = dict()
        for item in data:
            d_dict[pd.to_datetime(item[0])] = float(item[1].replace(',',''))
        qt = pd.Series(d_dict)
        
        
        flow_units = '$m^{3}s^{-1}$'
        #fix_units() # not needed
        #qt = float(scale_factor) * qt
        qt = qt.resample('D', label='right', loffset='-12h').mean()
        
        # NEW to deal with the problem that when you request date from before
        # 18 months ago it gives the most recent data instead.
        
        dt0_actual = qt.index[0]
        dt0_requested = days[0]
        import numpy as np
        if np.abs((dt0_actual - dt0_requested).days) >= 1:
            memo = 'That date range was not available'
            qt = ''
        else:
            got_data = True
            memo = 'success'
            
    except:
        memo = 'problem parsing data from soup'
    memo = (memo + ' EC')
    return qt, memo # NEW

def get_ec_data_historical(year):
    # gets Environment Canada data, using code cribbed from:
    #https://bitbucket.org/douglatornell/ecget/src/
    # NOTE: this will get data up through the end of 2015.
    import requests
    import bs4
    try:
        params = {
            'mode': 'Table',
            'type': 'h2oArc',
            'stn': '08MF005',
            'dataType': 'Daily',
            'parameterType': 'Flow',
            'year': str(year),
            'y1Max': '1',
            'y1Min': '1',
        }
        DATA_URL = 'http://wateroffice.ec.gc.ca/report/historical_e.html'
        DISCLAIMER_COOKIE = {'disclaimer': 'agree'}
        response = requests.get(DATA_URL, params=params,
                                cookies=DISCLAIMER_COOKIE)
        soup = bs4.BeautifulSoup(response.content, 'lxml')
        if (str(year) + ' Daily Discharge') in soup.text:
            # we do the test above because the website will return the
            # table for the most recent year if the requested year is
            # missing
            table = soup.find('table')
            table_body = table.find('tbody')
            rows = table_body.find_all('tr')
            d_dict = dict()
            for row in rows:
                # what day is it
                for ele in row.find_all('th'):
                    day = int(ele.text.strip())
                cols = row.find_all('td')
                cols = [ele.text.strip() for ele in cols] # a list of strings
                imo = 1
                for item in cols:
                    if len(item) == 0:
                        pass
                    else:
                        this_data = float(item.split()[0].replace(',',''))
                        # the split call is to remove trailing letters
                        # that occasionally appear after the data
                        d_dict[datetime(year, imo, day, 12, 0, 0)] = this_data
                    imo += 1
            qt = pd.Series(d_dict)
            flow_units = '$m^{3}s^{-1}$'
            #fix_units() # not needed
            #qt = float(scale_factor) * qt
            got_data = True
            memo = 'success'
        else:
            memo = 'That year was not available'
            qt = '' # NEW
    except:
        memo = 'problem parsing data from soup'
    memo = (memo + ' EC')
    return qt, memo # NEW

# driver
for year in range(2015, 2019):
    
    days = (datetime(year,1,1), datetime(year,1,5))
    print('\nYear = ' + str(year))
    
    print('REALTIME')
    qt, memo = get_ec_data(days)
    if len(qt) > 0:
        print(' ' + memo + ' ' + str(qt.index[0])+ ' ' + str(qt.index[-1]))
    else:
        print(memo)
        
    print('HISTORICAL')
    qth, memoh = get_ec_data_historical(year)
    if len(qth) > 0:
        print(' ' + memoh + ' ' + str(qth.index[0])+ ' ' + str(qth.index[-1]))
    else:
        print(memoh)
        
    sys.stdout.flush()


