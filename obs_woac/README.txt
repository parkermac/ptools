Programs for processing WOAC cruise data.  These are CTD+Bottle casts of data provided by Simone Alin (11/2018).  The code is structured to be much like that in obs_ecy and obs_collias.

======================================================================
* process_data.py

Input: ptools_data/woac/raw/WOAC_data_9-7-2018_dataToParker.xlsx

Output: ptools_data/woac/sta_df.p
		ptools_data/woac/Casts_2017.p
		
sta_df is a DataFrame with info like:
        Station        Cruise Latitude Longitude             Datetime
castnum                                                              
0.0          28       CAB1065  47.7108  -122.454  2017-04-06 17:27:39
1.0         28b       CAB1065  47.7042  -122.452  2017-04-06 18:18:01
2.0           5       CAB1065   47.885  -122.373  2017-04-06 21:24:59
3.0           7       CAB1065  47.9847  -122.623  2017-04-06 23:18:31
4.0           8       CAB1065  47.8989   -122.61  2017-04-07 00:24:22
5.0          12       CAB1065  47.4266  -123.108  2017-04-07 19:29:19
6.0         402       CAB1065  47.3571  -123.021  2017-04-07 21:38:10
7.0          11       CAB1065   47.373  -123.131  2017-04-07 22:30:12
8.0          13       CAB1065  47.5494  -123.004  2017-04-07 00:05:10

Note that the castnum index is a much cleaner way of identifying casts, and this carries over to the data in Casts_2017.

Note that when I extract casts from the model I append 'WOAC' to the station name (LiveOcean/x_cast).

Casts_2017 is a DataFrame with data like:
      Cruise  Longitude  Latitude Station  Pressure (dbar)   ...     Omega Ar     DO (uM)       Z (m)            Datetime  castnum
0    CAB1065 -122.45541  47.70757      28          183.761   ...     0.833404  252.984516 -182.142790 2017-04-06 17:20:07      0.0
1    CAB1065 -122.45550  47.70774      28          184.378   ...          NaN  252.806723 -182.754081 2017-04-06 17:20:47      0.0
2    CAB1065 -122.45471  47.71024      28          144.985   ...     0.849518  258.225945 -143.721699 2017-04-06 17:25:41      0.0
3    CAB1065 -122.45382  47.71103      28          117.492   ...          NaN  261.564385 -116.475996 2017-04-06 17:27:37      0.0
4    CAB1065 -122.45382  47.71103      28          117.598   ...          NaN  261.545452 -116.581049 2017-04-06 17:27:38      0.0
5    CAB1065 -122.45382  47.71104      28          117.562   ...          NaN  261.480525 -116.545371 2017-04-06 17:27:39      0.0
6    CAB1065 -122.45382  47.71104      28          117.514   ...          NaN  261.446276 -116.497799 2017-04-06 17:27:39      0.0
7    CAB1065 -122.45422  47.71141      28          113.366   ...     0.910803  261.921844 -112.386794 2017-04-06 17:30:06      0.0
8    CAB1065 -122.45422  47.71197      28           80.693   ...     0.873082  263.524334  -80.002335 2017-04-06 17:31:39      0.0
9    CAB1065 -122.45408  47.71254      28           51.645   ...     0.779552  262.735246  -51.206563 2017-04-06 17:33:14      0.0
10   CAB1065 -122.45402  47.71294      28           29.944   ...     0.779128  265.136759  -29.691353 2017-04-06 17:34:26      0.0
11   CAB1065 -122.45306  47.70394     28b           20.972   ...     0.795994  266.878131  -20.795522 2017-04-06 18:15:33      1.0
12   CAB1065 -122.45266  47.70409     28b           10.321   ...     0.882316  291.453175  -10.234414 2017-04-06 18:16:55      1.0
and the full list of columns is:
['Cruise', 'Longitude', 'Latitude', 'Station', 'Pressure (dbar)',
       'Temp. (deg C)', 'Salinity', 'Sigma (kg m-3)', 'DO Flag', 'NO3 (uM)',
       'NO2 (uM)', 'NH4 (uM)', 'Chl (ug L-1)', 'TA1 (umol kg-1)',
       'DIC1 (umol kg-1)', 'TA1 Flag', 'DIC1 Flag', 'TA2 (umol kg-1)',
       'DIC2 (umol kg-1)', 'TA2 Flag', 'DIC2 Flag', 'pH', 'pCO2 (uatm)',
       'CO2 (umol kg-1)', 'HCO3- (umol kg-1)', 'CO3-- (umol kg-1)', 'Omega Ca',
       'Omega Ar', 'DO (uM)', 'Z (m)', 'Datetime', 'castnum']

======================================================================
* station_map.py

Input: ptools_data/woac/sta_df.p

Output: ptools_output/woac/station_map.png

======================================================================
* plot_casts.py

Input:	ptools_data/woac/sta_df.p
		ptools_data/woac/Casts_2017.p
		LiveOcean_output/cast/cas6_v3_lo8b/WOAC[Station]_[YYYY.MM.DD].nc

Output: one plot for each station in ptools_output/woac/ e.g. CAB1065_1_2017.04.09.png
		meaning [Cruise]_[Station]_[YYYY.MM.DD].png


