#!/usr/bin/env python
"""
Testing script to check shaplet decomposition and plotting
"""

#Project ToDo:
#   interactive region selector
#   polar shapelets
#   fit for nmax
#   selecting the best beta_0

import sys,os
import pyfits as pf
import numpy as n
import pylab as p
from scipy import optimize

#shapelet functions
import img, decomp, shapelet

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_IMAGE')
    o.set_description(__doc__)
    o.add_option('-r', '--region', dest='region', default=None,
        help='Region of image to decompose into shapelets, (xmin,xmax,ymin,ymax), default: None')
    o.add_option('-N', '--noise_region', dest='nregion', default=None,
        help='Region of image to use to create a noise map, if set to None the entire image is used with an iterative process to clip out the tails, (xmin,xmax,ymin,ymax), default: None')
    o.add_option('-n', '--nmax', dest='nmax', default=5, type='int',
        help='Size of coefficient dimensions for minimization fit, default: 5')
    o.add_option('-R', '--range', dest='nrange', default=10, type='int',
        help='Range of N values test in a brute force fit after minimization, default: 10')
    o.add_option('-b', '--beta', dest='beta', default=None, type='float',
        help='Set an initial beta value, default: None, guess is made based on image size')
    o.add_option('--xtol', dest='xtol', default=0.0001, type='float',
        help='Relative error in parameters acceptable for convergence, default: 0.0001')
    o.add_option('--ftol', dest='ftol', default=0.0001, type='float',
        help='Relative error in chi^2 function acceptable for convergence, default: 0.0001')
    o.add_option('--maxiter', dest='maxiter', default=10, type='int',
        help='Maximum number of iterations to perform, default: 10')
    o.add_option('--frac', dest='frac', default=1., type='float',
        help='Fractional radius of image to fit the centroid within, default: 1, the entire image')
    opts, args = o.parse_args(sys.argv[1:])

    fn=args[0]
    hdulist=pf.open(fn)
    im=hdulist[0].data
    im=im[0,0]
    im0=im
    if not (opts.region is None):
        extent=map(int, opts.region.split(','))
        im=img.selPxRange(im,extent)

    #noise map
    if opts.nregion is None:
        nm=img.estimateNoiseMap(im)
    else:
        nextent=map(int, opts.nregion.split(','))
        nm=img.estimateNoiseMap(im0,region=nextent)
        nm=img.selPxRange(nm,extent)
    
    #select initial beta and xc
    if opts.beta==None:
        beta=decomp.initBeta(im)
    else:
        beta=[opts.beta,opts.beta]
    beta=decomp.initBeta2(im,frac=.2)
    #xc=img.centroid(im)
    xc=img.maxPos(im)
    nmax=[opts.nmax,opts.nmax]
    #centroid_frac=opts.frac
    #centroid_frac=.05
    print "beta0: (%f,%f)\tcentroid: (%f,%f)\tnmax: %i"%(beta[0],beta[1],xc[0],xc[1],nmax[0])

    #p.subplot(221)
    #p.imshow(im0)
    #if not (opts.region is None):
    #    extent=map(int, opts.region.split(','))
    #    rectx=[extent[0],extent[0],extent[1],extent[1],extent[0]]
    #    recty=[extent[2],extent[3],extent[3],extent[2],extent[2]]
    #    p.plot(rectx,recty,'k')
    #if not (opts.nregion is None):
    #    extent=map(int, opts.region.split(','))
    #    rectx=[nextent[0],nextent[0],nextent[1],nextent[1],nextent[0]]
    #    recty=[nextent[2],nextent[3],nextent[3],nextent[2],nextent[2]]
    #    p.plot(rectx,recty,'r')
    #p.colorbar()

    #p.subplot(222)
    #p.imshow(im)
    #p.text(xc[0],xc[1],'+')
    #p.colorbar()
    #
    #p.subplot(223)
    #p.imshow(nm)
    #p.colorbar()
    #
    #p.subplot(224)
    #p.hist(nm.flatten(),bins=100)
   
    #scipy optimize library downhill simplex minimization
    print 'Running minimization for beta and centroid...'
    xopt,fopt,iters,funcalls,warn,allvecs=optimize.fmin(decomp.chi2Func,[beta[0],beta[1],xc[0],xc[1]],args=(nmax,im,nm),xtol=opts.xtol,ftol=opts.ftol,maxiter=opts.maxiter,full_output=True,retall=True)
    #xopt,fopt,iters,funcalls,warn,allvecs=optimize.fmin(decomp.chi2betaFunc,[beta[0],beta[1]],args=(xc[0],xc[1],nmax,im,nm),xtol=opts.xtol,ftol=opts.ftol,maxiter=opts.maxiter,full_output=True,retall=True)
    print '\tDone'
    #scipy optimize brute force over a range of N values
    nrange=opts.nrange   #number of steps in N for the brute force
    if nmax[0]-int(nrange/2) > 0: n0=nmax[0]-int(nrange/2)
    else: n0=1
    n1=nrange+n0
    print 'Running brute force for size of N on range [%i:%i]...'%(n0,n1)
    x0=optimize.brute(decomp.chi2nmaxFunc,[n.s_[n0:n1:1]],args=(im,nm,[xopt[0],xopt[1]],[xopt[2],xopt[3]]),finish=None)
    #x0=optimize.brute(decomp.chi2nmaxFunc,[n.s_[n0:n1:1]],args=(im,nm,[xopt[0],xopt[1]],xc),finish=None)
    nmax=[int(x0),int(x0)]
    print '\tDone'

    #plot: data, model, residual: model-data, coeffs
    p.subplot(221)
    p.title('Image')
    p.imshow(im)
    p.colorbar()
    
    p.subplot(222)
    p.title('Model')
    rx=n.array(range(0,im.shape[0]),dtype=float)-xopt[2]
    ry=n.array(range(0,im.shape[1]),dtype=float)-xopt[3]
    #rx=n.array(range(0,im.shape[0]),dtype=float)-xc[0]
    #ry=n.array(range(0,im.shape[1]),dtype=float)-xc[1]
    bvals=decomp.genBasisMatrix([xopt[0],xopt[1]],nmax,rx,ry)
    coeffs=decomp.solveCoeffs(bvals,im)
    mdl=img.constructHermiteModel(bvals,coeffs,[xopt[2],xopt[3]],im.shape)
    #mdl=img.constructHermiteModel(bvals,coeffs,[xc[0],xc[1]],im.shape)
    p.imshow(mdl)
    p.text(xopt[2],xopt[3],'+')
    #p.text(xc[0],xc[1],'+')
    p.colorbar()
    
    p.subplot(223)
    p.title('Residual')
    res=im-mdl
    p.imshow(res)
    p.colorbar()

    p.subplot(224)
    p.title('Coefficents')
    coeffs=n.reshape(coeffs,nmax)
    p.pcolor(coeffs)
    p.colorbar()
    
    p.show()

