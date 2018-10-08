"""
Code to do initial processing of mooring data from Willapa Bay.
"""

import pandas as pd
import numpy as np
import pickle
from datetime import datetime, timedelta

# where the data is, and where results will be stored
dir0 = '../../ptools_data/willapa/'

fn = dir0 + '2017Bay Center pCO2AverageBurke.xlsx'

df = pd.read_excel(fn)

# the original file timestamp does not have the year, and so
# pandas reads it as 1900. Here we try to make the correct time
dt = datetime(2017,1,1) - datetime(1900,1,1)
df['Date'] = df['Julian Day'] + dt

df = df.set_index('Date')

df.pop('Julian Day')

df = df.rename(columns={'Ω Arag':'Omega'})

