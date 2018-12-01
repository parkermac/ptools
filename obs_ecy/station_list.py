#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Makes a text list of station info suitable for use in JavaScript.
"""

import pandas as pd

dir0 = '../../ptools_data/ecology/'
# load processed station info and data
sta_df = pd.read_pickle(dir0 + 'sta_df.p')

# add Canadian data
dir1 = '../../ptools_data/canada/'
# load processed station info and data
sta_df_ca = pd.read_pickle(dir1 + 'sta_df.p')

sta_df = pd.concat((sta_df, sta_df_ca), sort=False)

print('sta_list = [')
for station in sta_df.index:
    lon = sta_df.loc[station, 'Longitude']
    lat = sta_df.loc[station, 'Latitude']
    print('    {sta:"%s", lat:%0.4f, lng:%0.4f},' % (station, lat, lon))
print('];')
