#!/usr/bin/env python
"""
Testing script to check shaplet decomposition and plotting
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
    o.add_option('-n', '--nmax', dest='nmax', default='5',
        help='Size of coefficient dimensions for minimization fit, can be two values i.e. \'4,5\', default: 5')
    o.add_option('-B', '--brute', dest='brute', default=10, type='int',
        help='Maximum basis function order to use when running brute force method, default: 10')
    o.add_option('-b', '--beta', dest='beta', default=None, type='float',
        help='Set an initial beta value, default: None, guess is made based on image size')
    o.add_option('-o', '--outfile', dest='ofn', default='tempCart.coeff',
        help='Coefficients output filename, default: tempCart.coeff')
    o.add_option('--xtol', dest='xtol', default=0.0001, type='float',
        help='Relative error in parameters acceptable for convergence, default: 0.0001')
    o.add_option('--ftol', dest='ftol', default=0.0001, type='float',
        help='Relative error in chi^2 function acceptable for convergence, default: 0.0001')
    o.add_option('--maxiter', dest='maxiter', default=10, type='int',
        help='Maximum number of iterations to perform, default: 10')
    o.add_option('--frac', dest='frac', default=1., type='float',
        help='Fractional radius of image to fit the centroid within, default: 1, the entire image')
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
    
    nmax=opts.nmax.split(',')
    if len(nmax)==1:
        nmax=[int(nmax[0])+1,int(nmax[0])+1]
    else:
        nmax=[int(nmax[0])+1,int(nmax[1])+1]
    print "beta0: (%f,%f)\tcentroid: (%f,%f)\tnmax: %i"%(beta[0],beta[1],xc[0],xc[1],nmax[0])

    #scipy optimize library downhill simplex minimization
    print 'Running minimization for beta and centroid...'
    xopt,fopt,iters,funcalls,warn,allvecs=optimize.fmin(shapelets.decomp.chi2Func,[beta[0],beta[1],xc[0],xc[1]],args=(nmax,im,nm),xtol=opts.xtol,ftol=opts.ftol,maxiter=opts.maxiter,full_output=True,retall=True)
    print '\tDone'
    
    beta0=[xopt[0],xopt[1]]
    xc0=[xopt[2],xopt[3]]

    #scipy optimize brute force over a range of N values
    n0=1
    n1=opts.brute+1
    print 'Running brute force for size of N on range [%i:%i]...'%(n0,n1-1)
    x0=optimize.brute(shapelets.decomp.chi2nmaxFunc,[n.s_[n0:n1:1]],args=(im,nm,[xopt[0],xopt[1]],[xopt[2],xopt[3]]),finish=None)
    nmax0=[int(x0),int(x0)]
    print '\tDone'

    #plot: data, model, residual: model-data, coeffs
    p.subplot(221)
    p.title('Image')
    p.imshow(im)
    p.colorbar()
    
    p.subplot(222)
    p.title('Model')
    rx=n.array(range(0,im.shape[0]),dtype=float)-xc0[0]
    ry=n.array(range(0,im.shape[1]),dtype=float)-xc0[1]
    bvals=shapelets.decomp.genBasisMatrix(beta0,nmax0,rx,ry)
    coeffs=shapelets.decomp.solveCoeffs(bvals,im)
    mdl=shapelets.img.constructModel(bvals,coeffs,xc0,im.shape)
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
    sqCoeffs=n.reshape(coeffs,nmax0)
    p.pcolor(sqCoeffs)
    p.colorbar()
    
    ofn=opts.ofn
    print 'Writing to file:',ofn
    shapelets.fileio.writeHermiteCoeffs(ofn,coeffs,xc0,im.shape,beta0,nmax0,info=ifn)
    
    p.show()

