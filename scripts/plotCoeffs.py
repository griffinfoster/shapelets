#!/usr/bin/env python
"""
Plot a Hermite shapelet coefficient file
"""

import sys,os
import numpy as n
import pylab as p
import shapelets

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] COEFF_FILE')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    d=shapelets.fileio.readCoeffs(ifn)

    if d['mode'].startswith('herm'):
        rx=n.array(range(0,d['size'][0]),dtype=float)-d['xc'][0]
        ry=n.array(range(0,d['size'][1]),dtype=float)-d['xc'][1]
        bvals=shapelets.decomp.genBasisMatrix(d['beta'],d['norder'],rx,ry)
        mdl=shapelets.img.constructModel(bvals,d['coeffs'],d['xc'],d['size'])
        p.suptitle('Hermite')
    elif d['mode'].startswith('lag'):
        r0,th0=shapelets.shapelet.polarArray(d['xc'],d['size'])
        bvals=shapelets.decomp.genPolarBasisMatrix(d['beta'],d['norder'],r0,th0)
        mdl=n.abs(shapelets.img.constructModel(bvals,d['coeffs'],d['xc'],d['size']))
        coeffs=d['coeffs']
        cimR=shapelets.img.polarCoeffImg(coeffs.real,d['norder'])
        cimI=shapelets.img.polarCoeffImg(coeffs.imag,d['norder'])
        cimI=n.fliplr(cimI)
        cim=n.concatenate((cimR,cimI),axis=1)
        p.suptitle('Laguerre')
    
    p.subplot(121)
    p.title('Model')
    p.imshow(mdl)
    p.text(d['xc'][1],d['xc'][0],'+')
    p.colorbar()

    p.subplot(122)
    p.title('Coefficents')
    if d['mode'].startswith('herm'):
        coeffs=n.reshape(d['coeffs'],d['norder'])
        p.pcolor(coeffs)
    elif d['mode'].startswith('lag'):
        p.pcolor(cim)
    p.colorbar()
    
    p.show()
