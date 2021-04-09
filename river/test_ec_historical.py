"""
Code to test pulling Environment Canada historical river data from tables.

Note: in the new loenv environment I had to add bs4 and lxml to the loenv.yml
to get this to work.  I also added the test for '-' in the data, which they use
for days that don't exist.

RESULT: On 2021.03.29 this returned data through the end of 2019.
"""

import requests
import bs4
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt

plt.close('all')
for year in range(2015, 2021):
    print(str(year).center(60,'-'))

    # gets Environment Canada data, using code cribbed from:
    #https://bitbucket.org/douglatornell/ecget/src/
    # NOTE: this will get data up through the end of 2015.

    params = {
        'mode': 'Table',
        'type': 'h2oArc',
        'stn': '08MF005', # Fraser River
        'dataType': 'Daily',
        'parameterType': 'Flow',
        'year': str(year),
        'y1Max': '1',
        'y1Min': '1',
    }
    DATA_URL = 'http://wateroffice.ec.gc.ca/report/historical_e.html'
    DISCLAIMER_COOKIE = {'disclaimer': 'agree'}
    response = requests.get(DATA_URL, params=params, cookies=DISCLAIMER_COOKIE)
    soup = bs4.BeautifulSoup(response.content, 'lxml')
    # Try print(soup.prettify()) to have a look.  Works for any of the elements as well.
    if (str(year) + ' Daily Discharge') in soup.text:
        # we do the test above because the website will return the
        # table for the most recent year if the requested year is
        # missing
        table = soup.find('table')
        table_body = table.find('tbody')
        # table_body is just the data part of the table, with "headers" for the day of month
        # at the start of each row, and the columns are the months (no headers)
        rows = table_body.find_all('tr')
        d_dict = dict()
        # the table_body is arranged with rows = day of month, and column = month
        for row in rows:
            # what day is it
            for ele in row.find_all('th'):
                day = int(ele.text.strip())
            cols = row.find_all('td')
            cols = [ele.text.strip() for ele in cols] # a list of strings
            imo = 1
            for item in cols:
                if (len(item) == 0) or (item=='-'):
                    pass
                else:
                    this_data = float(item.split()[0].replace(',',''))
                    # the split call is to remove trailing letters
                    # that occasionally appear after the data
                    d_dict[datetime(year, imo, day, 12, 0, 0)] = this_data
                imo += 1
        qt = pd.Series(d_dict).sort_index()
        print(' + got data')
        qt.plot()
    else:
        print(' - no data')
plt.show()
