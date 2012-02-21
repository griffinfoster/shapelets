#!/usr/bin/env python
"""
Return complex correlation data from a set of shapelet coefficients and UV positions
"""

import sys,os
import pyfits as pf
import numpy as n
import pylab as p
from scipy import optimize
import shapelets

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] COEFF_FILE')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    #load cofficients file
    ifn=args[0]
    d=shapelets.fileio.readCoeffs(ifn)

    #generate inverse basis functions
    print 'generating basis functions...',
    sys.stdout.flush()
    if d['mode'].startswith('herm'):
        bfs=shapelets.shapelet.ftHermiteBasis(d['beta'],d['norder'])
    elif d['mode'].startswith('lag'):
        bfs=ftLaguerreBasis(d['beta'],d['norder'])
    print 'done'
    
    #generate random UV samples
    print 'generating random UV samples...',
    sys.stdout.flush()
    nsamples=(100,4)
    scalef=[100.,70.]
    uu=n.random.random(nsamples)*scalef[0]
    vv=n.random.random(nsamples)*scalef[1]
    print 'done'
    
    #compute the UV complex correlation from coefficients and basis functions
    print 'generating UV correlations...',
    sys.stdout.flush()
    if d['mode'].startswith('herm'):
        corr=shapelets.uv.computeHermiteUV(bfs,d['coeffs'],uu,vv)
    elif d['mode'].startswith('lag'):
        corr=shapelets.uv.computeLaguereUV(bfs,d['coeffs'],uu,vv)
    print 'done'

    #return an array of complex amplitudes
    p.subplot(121)
    p.scatter(uu,vv,c=corr.real,edgecolors=None)
    p.colorbar()
    p.subplot(122)
    p.scatter(uu,vv,c=corr.imag,edgecolors=None)
    p.colorbar()

    p.show()

