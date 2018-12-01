This code gets river data from USGS and Environment Canada sources.  It focuses on long records that are relevant to the SSMSP work.

===================================================================
* make_long_historical.py extracts maximum length flow records for a few selected USGS river.

Input: Web data from USGS, ssmsp/river_class.py & ssmsp/river_info.csv (for gauge numbers and ungauged flow factors)

Output: ptools_output/river_long/[river name].p

===================================================================
* plot_long_historical.py plots the records created by make_long_historical.py.

Input: ptools_output/river_long/[river name].p

Output: ptools_output/river_long/Historical_[river name].png
