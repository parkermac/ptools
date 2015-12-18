import matplotlib.pyplot as plt 
plt.close()
fig1 = plt.figure()

ax1 = fig1.add_subplot(1, 1, 1)

y = []
yy = []
x = range(-15, 16)
for z in x:
    y.append(z*z + 2) 
    yy.append(.5*z*z + 2)
    
ax1.plot(x,y,'-ob')
ax1.plot(x,yy,'-og')

ax1.grid()
plt.show()
