#!/usr/bin/env python
"""
Plot a shapelet coefficient file
"""

import sys,os
import numpy as n
import pylab as p

#shapelet functions
import img, decomp, shapelet, fileio

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] COEFF_FILE')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    d=fileio.readHermiteCoeffs(ifn)

    rx=n.array(range(0,d['size'][0]),dtype=float)-d['xc'][0]
    ry=n.array(range(0,d['size'][1]),dtype=float)-d['xc'][1]
    bvals=decomp.genBasisMatrix(d['beta'],d['norder'],rx,ry)
    mdl=img.constructHermiteModel(bvals,d['coeffs'],d['xc'],d['size'])
    
    p.subplot(121)
    p.title('Model')
    p.imshow(mdl)
    p.text(d['xc'][0],d['xc'][1],'+')
    p.colorbar()

    p.subplot(122)
    p.title('Coefficents')
    coeffs=n.reshape(d['coeffs'],d['norder'])
    p.pcolor(coeffs)
    p.colorbar()
    
    p.show()
