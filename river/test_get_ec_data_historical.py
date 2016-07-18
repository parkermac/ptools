# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 08:56:56 2016

@author: PM5

Test of methods (used eventually in river_class.py) to get historical EC
river data

https://wateroffice.ec.gc.ca/report/report_e.html?mode=Table&type=h2oArc
&stn=08MF005&dataType=Daily&parameterType=Flow&year=2003&y1Max=1&y1Min=1
"""


import requests
import bs4
from datetime import datetime
import pandas as pd

year = 1995

ec_code = '08GB013'

params = {
    'mode': 'Table',
    'type': 'h2oArc',
    'stn': ec_code,
    'dataType': 'Daily',
    'parameterType': 'Flow',
    'year': str(year),
    'y1Max': '1',
    'y1Min': '1',
}
DATA_URL = 'http://wateroffice.ec.gc.ca/report/report_e.html'
DISCLAIMER_COOKIE = {'disclaimer': 'agree'}
response = requests.get(DATA_URL, params=params, cookies=DISCLAIMER_COOKIE)
soup = bs4.BeautifulSoup(response.content, 'lxml')
if (str(year) + ' Daily Discharge') in soup.text:
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
                this_data = float(item.split()[0]) # remove trailing letters
                #print(datetime(year, imo, day))
                #print(this_data)
                d_dict[datetime(year, imo, day, 12, 0, 0)] = this_data
            imo += 1
    qt = pd.Series(d_dict)
else:
    print('That year is not available')
