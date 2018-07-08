# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 15:28:16 2016

@author: PM5
"""

from bokeh.plotting import figure, output_file, show

import numpy as np

# prepare some data
x = np.arange(100)
y = x**2

# output to static HTML file
output_file("../../ptools_output/scratch/lines.html")

# create a new plot with a title and axis labels
p = figure(title="simple line example", x_axis_label='x', y_axis_label='y')

# add a line renderer with legend and line thickness
p.line(x, y, legend="Temp.", line_width=5)

# show the results
show(p)