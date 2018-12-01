This is for processing the Collias dataset of historical CTD-Bottle data for Puget Sound.  This was minly done as part of the SSMSP work.  It included programs that combine with modern Ecology data to create longer time series.

======================================================================
* process_station_info.py makes a DataFrame of station locations, etc.

Input: ptools_data/collias/raw/Historic_Stations.xlsx

Output: ptools_data/collias/sta_df.p, a pickled DataFrame, with data like:
                           Descrip   Latitude   Longitude
Station                                                  
ADM201                MIDDLE POINT  48.176667 -122.776667
ADM202               PORT TOWNSEND  48.136667 -122.686667

======================================================================
* station_map.py plots the station locations from sta_df

Input: ptools_data/collias/sta_df.p

Output: ptools_output/collias/station_map.png

======================================================================
* process_data.py is the initial pass through the data.  It renames variables, keeping only selected ones

Input: ptools_data/collias/raw/Historic_Results.xlsx

Output: ptools_data/collias/Bottles_[1932-1975].p, with data like:
     Station  Time (UTC)  Temp. (deg C)  Salinity  DO (mg L-1)  ...   Letters  Numbers  Num0       Date Z (m)
0     ADM201      1849.0           7.76     31.16          NaN  ...       ADM      201     2 1952-04-17   -45
1     ADM201      1849.0           7.85     30.84          NaN  ...       ADM      201     2 1952-04-17   -30
2     ADM201      1849.0           7.96     30.52          NaN  ...       ADM      201     2 1952-04-17   -20
3     ADM201      1849.0           8.18     30.00          NaN  ...       ADM      201     2 1952-04-17   -10
4     ADM201      1849.0           8.19     29.99          NaN  ...       ADM      201     2 1952-04-17    -5
5     ADM201      1849.0           8.18     29.98         8.59  ...       ADM      201     2 1952-04-17     0

This is listing a single cast, but all the casts for the year are there.

Station is split into the columns Letters and Numbers, and the column Num0 is the first digit of Numbers, which designates the "region" we are in:
Regions (the first digit is stored ar Num0 in Bottles):
        100        Strait of Juan de Fuca
        200        Admiralty Inlet
        300        Puget Sound Basin
        400        Southern Puget Sound
        500        Hood Canal
        600        Whidbey Basin
        700        North Sound
        800        San Juan Island Passages
		
The full list of columns in this DataFrame is:
['Station', 'Time (UTC)', 'Temp. (deg C)', 'Salinity', 'DO (mg L-1)',
       'NO3 (mg L-1)', 'NO2 (mg L-1)', 'SiOH4 (mg L-1)', 'Letters', 'Numbers',
       'Num0', 'Date', 'Z (m)']
	   
I make "Z (m)" using the field "depth5" in the original, which I got from Skip Albertson.  I assume this is the actual recorded depth to the nearest 5 meters.  Given the interpolation I do later I could probably just use the actual depth, but I doubt it would make any difference.

======================================================================
* bin_by_region.py is the next step in the processing.  It gathers all the casts in a Region, for all years, and then for each cast interpolates in the vertical to just have data at [-30, -10, 0] m.

Input: Bottles_[1932-1975].p from above

Output: ptools_output/collias/region_[1-8].py which are DataFrames containing data like:
            Temp. (deg C)  Salinity  DO (mg L-1)  Z (m)  NO3 (uM)
Date                                                             
1952-08-04            NaN       NaN          NaN    -30       NaN
1952-08-04          12.73     29.93         9.95    -10       NaN
1952-08-04          12.49     29.56         8.61      0       NaN
1952-08-05            NaN       NaN          NaN    -30       NaN
1952-08-05          13.39     29.56          NaN    -10       NaN
1952-08-05          13.18     29.10          NaN      0       NaN

Issue: I have dropped all station information, just keeping all casts (just 3 depths) in a given Region.  The index is by Date.  Is there any possible problem with having multiple casts on the same date? My code uses "set(date_list)" and I assume this represents a SINGLE cast, but isn't it likely that they made more than one cast a day in a Region!? Answer: I think this is not a problem becasue the code is first pulling out all the data for a single station, and then identifying casts by assuming they have unique dates.  This is fine.

======================================================================
* bin_by_month_depth.py takes the data in each Region and forms sub-DataFrames of time series that average all the data in a region, in a month, in one of the three depth ranges.

Input: ptools_output/collias/region_[1-8].py &
	   ptools_output/ecology/region_[1-8].py

Output: ptools_output/collias/region_month_z_means.py which is a dict: A_dict[region] = (A0, A10, A30), numbered by depth, and each with data like (access each DataFrame using a command like A10 = A_dict[3][1])
           Temp. (deg C) Salinity DO (mg L-1) NO3 (mg L-1)
Date                                                      
1932-01-15           NaN      NaN         NaN          NaN
1932-02-15           NaN      NaN         NaN          NaN
1932-03-15           NaN      NaN         NaN          NaN
1932-04-15           NaN      NaN         NaN          NaN
1932-05-15           NaN      NaN         NaN          NaN
1932-06-15       13.3271  24.8157     8.97714          NaN

So this is a time series with 12 values per year, at mid month, over the years 1932-2017, which is 86 years! (although not all have data)

======================================================================
* plot_binned_region_series.py plots time series over the full record of the combined Collias and Ecology datasets.

Input: region_month_z_means.p (made by bin_by_month_depth.py)

Output: plots in ptools_output/collias/Series_[region name].png

======================================================================
* plot_binned_region_monthly.py plots monthly climatologies for three time periods from the combined Collias and Ecology datasets.

Input: region_month_z_means.p (made by bin_by_month_depth.py)

Output: plots in ptools_output/collias/Monthly_[region name].png

======================================================================
* plot_ts_vs_z.py is designed to explore the apparent changes over time that we saw in the plots from plot_binned_region_[series, monthly].py above (e.g. the surface of Hood Canal is fresher and warmer).  We just look at Temperature and Salinity because these are presumably the most robust signals.

Input: ptools_data/collias/Bottles_[1932-1975].p, and
	   ptools_data/ecology/Casts_[1999-2017].p
	   
Output: plots of T or s vs. z, showing both the raw data for the two datasets (Collias and Ecology) and binned by depth profiles in selected time periods.  These are saved as plots in ptools_output/collias/TSZ_[region name].png



