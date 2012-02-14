#!/usr/bin/env python
"""
Testing script to check polar shaplet decomposition and plotting
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
    o.set_usage('%prog [options] FITS_IMAGE')
    o.set_description(__doc__)
    o.add_option('-r', '--region', dest='region', default=None,
        help='Region of image to decompose into shapelets, (xmin,xmax,ymin,ymax), default: None')
    o.add_option('-N', '--noise_region', dest='nregion', default=None,
        help='Region of image to use to create a noise map, if set to None the entire image is used with an iterative process to clip out the tails, (xmin,xmax,ymin,ymax), default: None')
    o.add_option('-n', '--nmax', dest='nmax', default=7, type='int',
        help='Size of coefficient dimensions for minimization fit, default: 5')
    o.add_option('-B', '--brute', dest='brute', default=10, type='int',
        help='Maximum basis function order to use when running brute force method, default: 10')
    o.add_option('-b', '--beta', dest='beta', default=None, type='float',
        help='Set an initial beta value, default: None, guess is made based on image size')
    o.add_option('-o', '--outfile', dest='ofn', default='tempPolar.coeff',
        help='Coefficients output filename, default: tempPolar.coeff')
    o.add_option('--xtol', dest='xtol', default=0.0001, type='float',
        help='Relative error in parameters acceptable for convergence, default: 0.0001')
    o.add_option('--ftol', dest='ftol', default=0.0001, type='float',
        help='Relative error in chi^2 function acceptable for convergence, default: 0.0001')
    o.add_option('--maxiter', dest='maxiter', default=10, type='int',
        help='Maximum number of iterations to perform, default: 10')
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    im=shapelets.fileio.readFITS(args[0])
    im0=im
    if not (opts.region is None):
        extent=map(int, opts.region.split(','))
        im=shapelets.img.selPxRange(im,extent)

    #noise map
    if opts.nregion is None:
        nm=shapelets.img.estimateNoiseMap(im)
    else:
        nextent=map(int, opts.nregion.split(','))
        nm=shapelets.img.estimateNoiseMap(im0,region=nextent)
        nm=shapelets.img.selPxRange(nm,extent)
    
    #select initial beta and xc
    if opts.beta==None:
        beta=shapelets.decomp.initBeta(im)
    else:
        beta=[opts.beta,opts.beta]
    beta=shapelets.decomp.initBeta2(im,frac=.2)
    #xc=shapelets.img.centroid(im)
    xc=shapelets.img.maxPos(im)
    nmax=opts.nmax+1
    beta=(beta[0]+beta[1])/2.
    print "beta0: (%f)\tcentroid: (%f,%f)\tnmax: %i"%(beta,xc[0],xc[1],nmax)

    #scipy optimize library downhill simplex minimization
    print 'Running minimization for beta and centroid...'
    xopt,fopt,iters,funcalls,warn,allvecs=optimize.fmin(shapelets.decomp.chi2PolarFunc,[beta,xc[0],xc[1]],args=(nmax,im,nm),xtol=opts.xtol,ftol=opts.ftol,maxiter=opts.maxiter,full_output=True,retall=True)
    print '\tDone'

    xc0=[xopt[1],xopt[2]]
    beta0=xopt[0]
    
    #scipy optimize brute force over a range of N values
    n0=1
    n1=opts.brute+1
    print 'Running brute force for size of N on range [%i:%i]...'%(n0,n1-1)
    x0=optimize.brute(shapelets.decomp.chi2nmaxPolarFunc,[n.s_[n0:n1:1]],args=(im,nm,beta0,xc0),finish=None)
    print '\tDone'
    
    nmax0=int(x0)

    #plot: data, model, residual: model-data, coeffs
    p.subplot(221)
    p.title('Image')
    p.imshow(im)
    p.colorbar()
    
    p.subplot(222)
    p.title('Model')
    r0,th0=shapelets.shapelet.polarArray(xc0,im.shape)
    bvals=shapelets.decomp.genPolarBasisMatrix(beta0,nmax0,r0,th0)
    coeffs=shapelets.decomp.solveCoeffs(bvals,im)
    mdl=n.abs(shapelets.img.constructModel(bvals,coeffs,xc0,im.shape))
    p.imshow(mdl)
    p.text(xc0[1],xc0[0],'+')
    p.colorbar()
    
    p.subplot(223)
    p.title('Residual')
    res=im-mdl
    p.imshow(res)
    p.colorbar()
    
    p.subplot(224)
    p.title('Coefficents')
    cimR=shapelets.img.polarCoeffImg(coeffs.real,nmax0)
    cimI=shapelets.img.polarCoeffImg(coeffs.imag,nmax0)
    cimI=n.fliplr(cimI)
    cim=n.concatenate((cimR,cimI),axis=1)
    p.pcolor(cim)
    p.colorbar()

    ofn=opts.ofn
    print 'Writing to file:',ofn
    shapelets.fileio.writeLageurreCoeffs(ofn,coeffs,xc0,im.shape,beta0,nmax0,info=ifn)
    
    p.show()

