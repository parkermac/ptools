"""
Process station info from Collias casts, and
save the results for future use.

"""

import pandas as pd

# where the data is, and where results will be stored
dir0 = '../../ptools_data/collias/'

# file with station info (good for all years?)
sta_fn = dir0 + 'raw/Historic_Stations.xlsx'

# PROCESS STATION LOCATION INFO
sta_df = pd.read_excel(sta_fn)

col_dict = {'station':'Station', 'description':'Descrip',
    'latD':'Latitude', 'longD':'Longitude'}
    
sta_df = sta_df.rename(columns=col_dict)
sta_df = sta_df.set_index('Station')
sta_df = sta_df[['Descrip', 'Latitude', 'Longitude']]

# remove one station that has repeated, conflicting location info
sta_df = sta_df.drop(index='SKG628')

# save result
sta_df.to_pickle(dir0 + 'sta_df.p')

# gather unique basin identifiers
bas_list = []
for sta in sta_df.index:
    bas_list.append(sta[:3])
ubas_list = set(bas_list)
# RESULT:
# There are 605 stations, such as 'ADM201',
# and 52 3-letter idendifiers, such as 'ADM'.
