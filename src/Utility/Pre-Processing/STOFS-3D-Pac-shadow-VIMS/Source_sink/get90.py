import numpy as np
d = np.loadtxt('vsource.th.0')
p90 = np.percentile(d,90,axis=0)[:,np.newaxis].T # make into a row array
np.savetxt('vsource.th.p90',p90,fmt='%.0f')

