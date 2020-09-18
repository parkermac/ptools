"""
Code to get names of pandas line colors
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# make some lines to plot
x = np.linspace(0,10,300)

a = pd.DataFrame(index=x)

plt.close('all')
fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(111)

for ii in np.linspace(.1, 1,6):
    y = np.sin(x**ii)
    a.loc[:,ii] = y
    
a.plot(ax=ax)

for line in ax.get_lines():
    print(line.get_color())

plt.show()
