#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
sys.path.append('./shapelets/')
import matplotlib.pyplot as plt
import numpy as np
from time import time
import cshapelet as csh
import fshapelet as fsh
import shapelet as sh
import decomp
import cProfile, pstats
n = 16
m = 2
PB_DIM = (1000, 500)
xc = (800, 1300)
beta0 = 100.
#mask = (np.floor(np.random.rand(*PB_DIM)*1.1)).astype('bool').astype('uint8')
#print "Number of non null elements:", mask.sum(), 'out of', np.prod(PB_DIM)
mask = None
r, th = sh.polarArray(xc, PB_DIM)
print 'Problem dim:', PB_DIM
pr = cProfile.Profile()
pr.enable()
t = time()
#res0 = csh.polar_basis_L(r,n,m, beta0)*csh.polar_basis_Y(th, m,)
res0 = fsh.genPolarBasisMatrix(beta0, n, r, th)#{, mask=mask)
t = time() - t
pr.disable()
print 'Numexpr:',t


t = time()
#res1 = sh.polarDimBasis(n,m, beta0)(r,th)
res1 = decomp.genPolarBasisMatrix(beta0, n, r, th)
print 'Original version:', time() - t

ps = pstats.Stats(pr)
ps.sort_stats('cumulative').print_stats(20)

#mslice = 800 
#ax0 = plt.subplot(211)
#ax0.plot(np.abs(res0[mslice]), label='cython')
#ax0.plot(np.abs(res1[mslice]), label='python')
#ax0.legend()
##cs = plt.imshow(np.abs(res0/res1))
##plt.colorbar(cs)
#ax1 = plt.subplot(212)
#ax1.plot(np.abs(res0[mslice]/res1[mslice]), label='python')
#plt.show()

