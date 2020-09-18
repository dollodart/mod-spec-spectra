import numpy as np
mod_per = 300.
mod_freq = 1/mod_per
t = np.arange(0,5*mod_per,1) # time range, scan rate is 1 per second

cA = 1 + np.sin(mod_freq*t)
cAs = cA
cBs = cAs**2
cB = cBs

x = np.arange(1500,4000,1) # wavenumber range
sA = np.exp(-(1800-x)**2/10) + np.exp(-(3200-x)**2/40)
sAs = np.exp(-(2200-x)**2/10) + np.exp(-(2800-x)**2/40)
sBs = np.exp(-(2100-x)**2/10) + np.exp(-(2900-x)**2/40)
sB = np.exp(-(1600-x)**2/10) + np.exp(-(2600-x)**2/40)

c = np.array([cA,cAs,cBs,cB])
s = np.array([sA,sAs,sBs,sB]).transpose()
m = s @ c

#import matplotlib.pyplot as plt
for c,row in enumerate(m.transpose()):
    np.savetxt('example-data/data-t{}.csv'.format(c), 
            np.array([x,row]).transpose(),delimiter=',')
#    plt.plot(row,'k-')
#plt.show()
