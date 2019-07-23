"""
Code to iteratively solve for the flux fraction coefficients
at a channel segment.
"""

import numpy as np

# Specify the inflow (i) and outflow(o) transport and salinity.
# The index into each array corresponds to which section is it.

if True:
    # For this example we expect:
    #   Abest[0,0] = alpha34 = .4
    #   Abest[1,1] = alpha21 = .6
    Qi = np.array([25, 25])
    Qo = np.array([20, 30])
    Si = np.array([20,30])
    So = np.array([25, 25])
else:
    # Triple junction
    Qi = np.array([25, 25, 20])
    Qo = np.array([20, 30, 20])
    Si = np.array([20,30, 27])
    So = np.array([25, 25, 22])

# salt fluxes
Fi = Qi * Si
Fo = Qo * So
# normalize
qi = Qi / np.mean((Qi+Qo)/2)
qo = Qo / np.mean((Qi+Qo)/2)
fi = Fi / np.mean((Fi+Fo)/2)
fo = Fo / np.mean((Fi+Fo)/2)

ns = len(Qi) # number of sections

def get_err(ns, qi, qo, fi, fo, alo, ahi):
    # Create random transport fraction array.
    
    # The columns (index=1) tell which section this is
    # and the rows (index=0) tell which section the flow is
    # going to.  The sum of each column is one by construction,
    # ensuring that the net inflow is correct in each section,
    # but it does not ensure that the net outflow is correct.
    # This is done by the iterative solution.
    a = np.zeros((ns,ns))
    a[0,:] = alo + np.random.rand(ns)*(ahi-alo)
    if ns == 3:
        a[1,:] = np.random.rand(ns)*(1-a[0,:])
    a[-1,:] = 1 - np.sum(a[:-1,:], axis=0)
    # calculate the error
    qe = a.dot(qi) - qo
    fe = a.dot(fi) - fo
    return a, qe.std(), fe.std()
    
def run_trials(ns, nit, alo, ahi):
    A = np.ones((nit,ns,ns))
    qem = np.ones(nit)
    fem = np.ones(nit)
    for ii in range(nit):
        A[ii,:,:], qem[ii], fem[ii] = get_err(ns, qi, qo, fi, fo, alo, ahi)
    # find trials with the smallest combined error
    em = qem + fem
    iems = np.argsort(em)
    qems = qem[iems]
    fems = fem[iems]
    As = A[iems,:,:]
    return As, qems, fems

nit = 1000
ntop = int(nit/20)

# run a series of trials with iteratively better constrained guesses
for itr in range(5):
    if itr == 0:
        alo = np.zeros(ns)
        ahi = np.ones(ns)
    else:
        alo = np.min(As[:ntop,0,:],axis=0)
        ahi = np.max(As[:ntop,0,:],axis=0)
    print('-iteration = %d alo[0]=%0.2f ahi[0]=%0.2f' % (itr, alo[0], ahi[0]))
    As, qems, fems = run_trials(ns, nit, alo, ahi)
    

Abest = As[0,:,:]
print('Abest:')
print(Abest)
print('qems=%0.5f fems=%0.5f' % (qems[0], fems[0]))
