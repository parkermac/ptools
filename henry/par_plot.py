"""
Plots parabolas and othe polynomials.
"""
import numpy as np
import matplotlib.pyplot as plt

aa = 10
x = np.linspace(-aa, aa,1001)

y1 = x**2 - 1
y2 = 1 / y1
y3 = (x**3 + 7*x**2 - x + 1) / (x**2 + 2*x + 7)

plt.close()
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

ax.plot(x,y1,'-r',linewidth=2)
ax.plot(x,y2,'-b',linewidth=2)
ax.plot(x,y3,'-g',linewidth=2)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim((-aa, aa))
ax.set_ylim((-aa, aa))
ax.legend(['Line 1', 'Line 2', 'Line 3'])
ax.grid()

plt.show()