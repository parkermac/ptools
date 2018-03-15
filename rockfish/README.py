#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 16:18:09 2018

@author: pm7

README for the rockfish code.

PM_to_nc.py repackages all the daily pickle files into a single NetCDF file
for each experiment.

I ran it on fjord with pstep=10, meaning that all the files have about
10k particles instead of 100k (and are 1.7 GB instead of 17 GB).  Since
there are daily releases for the first 60 days of an experiment, that
means that there would be 10k/60 = 167 particles per release (or 1670 had
I kept the full file size).

PM_plot_nc.py is a handy tool for looking at the results of a single
experiment and developing processing and filtering tools.  It can only plot
about a hundred tracks at a time.

count_by_region.py is the real workhorse.  It figures out how many particles
from each experiment ended up in a region (defined by a polygon in 
ptools_data/polygons_rockfish) like Hood Canal.  It puts all the results
in a DataFrame.

NOTE: 2018.03.15 In response to this request from Kelly Andrews:
   I was thinking it would be helpful to know the variation around
   which basin each of the particles from each "born" date ended up in
   for each "age" between age 90 and 120 days. For example, looking at
   the spreadsheet I attached in the previous message, for each release
   day (1-60), I want to know how many of those 167 particles are in each
   of the basins at age 90-120. So, for each experiment, there should be
   1,860 rows of data (60 release days x 31 settling days).
   Does that make sense? Since we don't really know where in 90-120 larval
   duration period these guys will settle out, I think it's best to
   quantify a range or distribution across these days. And looking at the
   data across the release dates might give us some hint of whether there
   is variation in where the larvae end up depending on whether a female
   releases larvae early in the period or late.
   
   I made count_by_region_2.py.
                                                          
"""

