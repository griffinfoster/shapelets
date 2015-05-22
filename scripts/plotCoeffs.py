#!/usr/bin/env python
"""
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
import shapelets

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] COEFF_FILE')
    o.set_description(__doc__)
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save the figure, requires filename')
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    d=shapelets.fileio.readCoeffs(ifn)

    if d['mode'].startswith('herm'):
        #model
        rx=np.array(range(0,d['size'][0]),dtype=float)-d['xc'][0]
        ry=np.array(range(0,d['size'][1]),dtype=float)-d['xc'][1]
        xx,yy=shapelets.shapelet.xy2Grid(rx,ry)
        bvals=shapelets.decomp.genBasisMatrix(d['beta'],d['norder'],d['phi'],xx,yy)
        mdl=shapelets.img.constructModel(bvals,d['coeffs'],d['size'])

        #coeffs
        coeffs=np.reshape(d['coeffs'],d['norder'])

        plt.suptitle('Hermite')

    elif d['mode'].startswith('lag'):
        #model
        r0,th0=shapelets.shapelet.polarArray(d['xc'],d['size'])
        bvals=shapelets.decomp.genPolarBasisMatrix(d['beta'],d['norder'],d['phi'],r0,th0)
        mdl=np.abs(shapelets.img.constructModel(bvals,d['coeffs'],d['size']))

        #coeffs
        coeffs=d['coeffs']
        cimR=shapelets.img.polarCoeffImg(coeffs.real,d['norder'])
        cimI=shapelets.img.polarCoeffImg(coeffs.imag,d['norder'])
        cimI=np.fliplr(cimI)
        coeffs=np.concatenate((cimR,cimI),axis=1)

        plt.suptitle('Laguerre')
    
    plt.subplot(211)
    plt.title('Model')
    plt.imshow(mdl)
    plt.text(d['xc'][1],d['xc'][0],'+')
    plt.colorbar()

    plt.subplot(212)
    plt.title('Coefficents')
    plt.pcolor(coeffs)
    plt.colorbar()
    
    if not (opts.savefig is None):
        plt.savefig(opts.savefig)
    else: plt.show()
