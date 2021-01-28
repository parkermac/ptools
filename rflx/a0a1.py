"""
Code to explore approximate and exact versions of the efflux-reflux coefficients.

Following notes of 1/25/2021.

"""

s = 20
ds0 = 5
ds1 = 5

sout0 = s - ds0
sin0 = s
sout1 = s
sin1 = s + ds1

A0 = (sout0/sin0)*(sin1 - sin0)/(sin1 - sout0)
A1 = (sin1/sout1)*(sout1 - sout0)/(sin1 - sout0)

ds = (ds0 + ds1)/2

a0 = (ds1/ds - ds0*ds1/(s*ds))/2
a00 = (1 - ds/s)/2

a1 = (ds0/ds + ds0*ds1/(s*ds))/2
a11 = (1 + ds/s)/2

print('A0 = %0.2f, a0 = %0.2f, a00 = %0.2f' % (A0, a0, a00))
print('A1 = %0.2f, a1 = %0.2f, a11 = %0.2f' % (A1, a1, a11))
