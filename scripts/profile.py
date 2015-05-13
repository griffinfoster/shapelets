#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
sys.path.append('../src/')
import matplotlib.pyplot as plt
import numpy as np
from time import time
import cshapelet as csh
import fshapelet as fsh
import shapelet as sh
import decomp
import cProfile, pstats

USE_PROFILER = True

n = 16
m = 2
PB_DIM = (1000, 500)
xc = (800, 1300)
beta0 = 100.
mask = None
r, th = sh.polarArray(xc, PB_DIM)
print 'Problem dim:', PB_DIM

backend_list = {'original': decomp,
                'opt + numexpr': fsh,
                'opt + cython': csh
                 }

pr = {}

for backend, mod in backend_list.iteritems():
    if USE_PROFILER:
        pr[backend] = cProfile.Profile()
        pr[backend].enable()
    t = time()
    res0 = mod.genPolarBasisMatrix(beta0, n, r, th),
    t = time() - t
    if USE_PROFILER:
        pr[backend].disable()
    print backend,':', t


if USE_PROFILER:
    for backend in backend_list:
        print '='*80
        print ' '*30, backend
        print '='*80
        ps = pstats.Stats(pr[backend])
        ps.sort_stats('cumulative').print_stats(20)
