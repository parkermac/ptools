
# NOTE: this will get data up through today, but can only go back
# 18 months into the past.
# To get longer records use get_ec_data_historical.

import requests
import bs4
import pandas as pd

from datetime import datetime
days = (datetime(2017,1,1), datetime(2017,1,5))

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

