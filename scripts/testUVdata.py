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

    #generate image from coefficients
    print 'generating mdl image...',
    sys.stdout.flush()
    if d['mode'].startswith('herm'):
        rx=n.array(range(0,d['size'][0]),dtype=float)-d['xc'][0]
        ry=n.array(range(0,d['size'][1]),dtype=float)-d['xc'][1]
        bvals=shapelets.decomp.genBasisMatrix(d['beta'],d['norder'],rx,ry)
        mdl=shapelets.img.constructModel(bvals,d['coeffs'],d['xc'],d['size'])
    elif d['mode'].startswith('lag'):
        r0,th0=shapelets.shapelet.polarArray(d['xc'],d['size'])
        bvals=shapelets.decomp.genPolarBasisMatrix(d['beta'],d['norder'],r0,th0)
        mdl=n.abs(shapelets.img.constructModel(bvals,d['coeffs'],d['xc'],d['size']))
    print 'done'

    #generate random UV samples
    print 'generating random UV samples...',
    sys.stdout.flush()
    nsamples=(100,4)
    scalef=[100.,70.]
    uu=n.random.random(nsamples)*scalef[0]
    vv=n.random.random(nsamples)*scalef[1]
    print 'done'

    #cycle through samples and compute the DFT(inverse/scale factor?)
    print 'generating UV correlations...',
    sys.stdout.flush()
    rx=n.array(range(0,d['size'][0]),dtype=float)-d['xc'][0]
    ry=n.array(range(0,d['size'][1]),dtype=float)-d['xc'][1]
    rx0=n.reshape(n.tile(rx,len(ry)),(len(ry),len(rx)))
    ry0=n.reshape(n.tile(ry,len(rx)),(len(rx),len(ry)))
    corr=shapelets.dft.computeUV(mdl,rx0.T,ry0,uu,vv)
    print 'done'

    #return an array of complex amplitudes
    p.subplot(121)
    p.scatter(uu,vv,c=corr.real,edgecolors=None)
    p.colorbar()
    p.subplot(122)
    p.scatter(uu,vv,c=corr.imag,edgecolors=None)
    p.colorbar()

    p.show()

