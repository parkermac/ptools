
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D     # Import 3D plotting tools.
from matplotlib import cm
from scipy.special import jn                # Import Bessel function.

# Define grid of points.
points = np.linspace(-10, 10, 51)
X, Y = np.meshgrid(points, points)
R = np.sqrt(X**2 + Y**2)
Z = jn(0,R)

# Get lighting object for shading surface plots.
from matplotlib.colors import LightSource

# Get colormaps to use with lighting object.
from matplotlib import cm

# Create an instance of a LightSource and use it to illuminate the surface.
light = LightSource(90, 45)
illuminated_surface = light.shade(Z, cmap=cm.jet) # cm.coolwarm

rgb = np.ones((Z.shape[0], Z.shape[1], 3))
illuminated_surface = light.shade_rgb(rgb, Z)
ax = Axes3D(plt.figure())
# ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=False,
#                 facecolors=illuminated_surface)
                
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=False,
    facecolors=cm.jet(X/10))

                
plt.show()