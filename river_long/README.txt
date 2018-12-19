This code gets river data from USGS and Environment Canada sources.  It focuses on long records that are relevant to the SSMSP work.

The record length of these rivers is:
Skagit 1940-
Snohomish 1963-
Puyallup 1914-
Deschutes 1945-
Skokomish 1943- (bad after end of 2008)

===================================================================
* make_long_historical.py extracts maximum length flow records for a few selected USGS river.

Input: Web data from USGS, ssmsp/river_class.py & ssmsp/river_info.csv (for gauge numbers and ungauged flow factors)

Output: ptools_output/river_long/[river name].p

===================================================================
* plot_long_historical.py plots the records created by make_long_historical.py.

Input: ptools_output/river_long/[river name].p

Output: ptools_output/river_long/Historical_[river name].png
