#!/usr/bin/env python
"""
Plot a Shaplet coefficient file
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
import shapelets
import shapelets.phs

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] COEFF_FILE')
    o.set_description(__doc__)
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save the figure, requires filename')
    o.add_option('-f', '--ft', dest='fourier', action='store_true',
        help='Perform a Fourier transform on the basis functions')
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    d=shapelets.fileio.readCoeffs(ifn)

    ry=np.array(range(0,d['size'][0]),dtype=float)-d['xc'][0]
    rx=np.array(range(0,d['size'][1]),dtype=float)-d['xc'][1]

    print 'RA:', shapelets.phs.rad2hmsdms.rad2hmsdms(d['ra'],Type="ra",deg=True)
    print 'Dec:', shapelets.phs.rad2hmsdms.rad2hmsdms(d['dec'],Type="dec",deg=True)

    if d['mode'].startswith('herm'):
        #model
        yy,xx=shapelets.shapelet.xy2Grid(ry,rx)
        bvals=shapelets.decomp.genBasisMatrix(d['beta'],d['norder'],d['phi'],yy,xx,fourier=opts.fourier)
        mdl=shapelets.img.constructModel(bvals,d['coeffs'],d['size'])

        #coeffs
        coeffs=np.reshape(d['coeffs'],d['norder'])

        if opts.fourier:
            ptitle='Fourier(Hermite)'
        else:
            ptitle='Hermite'

    elif d['mode'].startswith('lag'):
        #model
        r0,th0=shapelets.shapelet.polarArray(d['xc'],d['size'])
        bvals=shapelets.decomp.genPolarBasisMatrix(d['beta'],d['norder'],d['phi'],r0,th0,fourier=opts.fourier)
        mdl=np.abs(shapelets.img.constructModel(bvals,d['coeffs'],d['size']))

        #coeffs
        coeffs=d['coeffs']
        cimR=shapelets.img.polarCoeffImg(coeffs.real,d['norder'])
        cimI=shapelets.img.polarCoeffImg(coeffs.imag,d['norder'])
        cimI=np.fliplr(cimI)
        coeffs=np.concatenate((cimR,cimI),axis=1)

        if opts.fourier:
            ptitle='Fourier(Laguerre)'
        else:
            ptitle='Laguerre'
    
    print shapelets.measure.all(d['coeffs'],d['beta'],d['norder'],mode=d['mode'])

    plt.suptitle(ptitle)
    plt.subplot(121)
    plt.title('Model')
    plt.imshow(np.abs(mdl),interpolation='nearest',extent=(rx[0],rx[-1],ry[-1],ry[0]))
    plt.text(0,0,'+',horizontalalignment='center',verticalalignment='center')
    plt.colorbar()

    plt.subplot(122)
    plt.title('Coefficients')
    #plt.pcolor(coeffs)
    plt.imshow(coeffs,interpolation='nearest',origin='lower')
    plt.colorbar()
    
    if not (opts.savefig is None):
        plt.savefig(opts.savefig)
    else: plt.show()

