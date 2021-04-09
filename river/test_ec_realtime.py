"""
Code to test the get_ec_data() function for current river data.
"""

import requests
import bs4
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

plt.close('all')
for year in range(2020,2021):
    print(str(year).center(60,'-'))
    days = (datetime(year,1,1), datetime(year,12,31))
    # gets Environment Canada data, using code originally cribbed from:
    #https://bitbucket.org/douglatornell/ecget/src/
    # NOTE: this will get data up through today, but can only go back
    # 18 months into the past.
    # To get longer records use get_ec_data_historical.
    try:
        PARAM_IDS = {'discharge': 47,}

        params = {
            'mode': 'Table',
            'type': 'realTime',
            'prm1': PARAM_IDS['discharge'],
            'prm2': -1,
            'stn': '08MF005', # Fraser River
            'startDate': days[0].strftime('%Y-%m-%d'),
            'endDate': days[1].strftime('%Y-%m-%d'),
        }
        DATA_URL = 'http://wateroffice.ec.gc.ca/report/real_time_e.html'
        DISCLAIMER_COOKIE = {'disclaimer': 'agree'}
        response = requests.get(DATA_URL, params=params, cookies=DISCLAIMER_COOKIE)
        soup = bs4.BeautifulSoup(response.content, 'lxml')
        table = soup.find('table')
        # If the data request was out of the time range (most recent 18 months, I think)
        # then "table" will be of type "NoneType" and the next command will
        # throw an exception.
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
        qth = pd.Series(d_dict)
        qt = qth.resample('D').mean() # the raw data is hourly
        qt.index += timedelta(days=0.5) # center data on noon
        print(' + got data')
        qt.plot()
    except Exception as e:
        print(' - no data')
        print(e)
        
plt.show()



