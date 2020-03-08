"""
Test of text placement and rotation.

RESULT: it is much easier to place and rotate text when using ha=va='center'.
In this case the ceanter of the text and the center of rotation are the same.

If you use 'left' and 'bottom' (the defaults) then the center of rotation is
on the lower right, leading to unexpected results.
"""

import matplotlib.pyplot as plt

plt.close('all')
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111)

ax_hand = ax.axis([-1,1,-1,1])

ii = 0
cvec = list('rbg')
for r in [0,-45, -90]:
    print(r)
    c = cvec[ii]
    line_hand = ax.plot(0,0,marker='*',color=c)
    if False:
        ha='left'
        va='bottom'
    else:
        ha = 'center'
        va = 'center'
    text_hand = ax.text(0,0,'Now is the time\nfor all good men',
        ha=ha, va=va,
        rotation=r, alpha=.5, size=50, color=c)
    ii+=1
    
plt.show()