"""
Make your own colormaps.

"""

# imports
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np

# make a data field to plot
x = np.linspace(-10,10,num=500)
y = x.copy()
X, Y = np.meshgrid(x,y)
R = np.sqrt(X**2 + Y**2)
Z = np.cos(3*R) * np.exp(-R/5)

# Here we make a colormap from a list of colors:
cmap = LinearSegmentedColormap.from_list('pm',['r','b',(0,0,0),'m'])
# The first argument 'pm' is the name assigned to the colormap object
# and it is returned by cmap.name.
# The second argument is a list of color names or (R,B,G) triple.
# By default the cmap object has 256 color steps, and you can see the (R,G,B,A)
# of each one by passing a number (0-255) as an argument:
# cmap(10) => (0.8823529411764706, 0.0, 0.11764705882352941, 1.0)

# You can use the colors modulr get more info on RGB values, e.g.
# colors.to_rgb('b') gives the RGB triple for blue: (0.0, 0.0, 1.0)

# Here we make a colormap by concatenating two colormaps:
top = cm.get_cmap('Blues', 128)
bottom = cm.get_cmap('rainbow', 128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
cmap2 = ListedColormap(newcolors, name='br')
# Note: this call top(np.linspace(0, 1, 128))
# returns a numpy array of RGBA with shape (128,4)

plt.close('all')
fig = plt.figure(figsize=(20,7))

ax = fig.add_subplot(121)
cs = ax.pcolormesh(X,Y,Z, cmap=cmap)
ax.axis('square')
fig.colorbar(cs)
ax.set_title('Colormap from list')

ax = fig.add_subplot(122)
cs = ax.pcolormesh(X,Y,Z, cmap=cmap2)
ax.axis('square')
fig.colorbar(cs)
ax.set_title('Colormap from colormaps')

plt.show()

