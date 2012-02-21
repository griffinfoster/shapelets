#!/usr/bin/env python
"""
Solve for shapelet coefficients based on beta, xc, and n_max
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
    o.add_option('-x', '--xc', dest='xc', default=None,
        help='set a x,y position for centroid, default: center of image/region')
    o.add_option('-m', '--mode', dest='mode', default='cart',
        help='Set the shapelet mode, cartesian or polar, default: cartesian')
    o.add_option('-n', '--nmax', dest='nmax', default='5',
        help='Size of coefficient dimensions, can be two values i.e. \'4,5\', default: 5')
    o.add_option('-b', '--beta', dest='beta', default=None,
        help='Beta value, can be two values i.e. \'25.0,30.5\', default: None, guess is made based on image size')
    o.add_option('-o', '--outfile', dest='ofn', default='temp.coeff',
        help='Coefficients output filename, default: temp.coeff')
    o.add_option('--max', dest='max_pos', action="store_true", default=False,
        help='Override centroid position to be the position of max intensity')
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    im=shapelets.fileio.readFITS(args[0])
    im0=im
    extent=[0,im.shape[0],0,im.shape[1]]
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
        beta=shapelets.decomp.initBeta2(im,frac=.2)
    else:
        beta=map(float,opts.beta.split(','))
        if len(beta)==1:
            beta=[beta[0],beta[0]]
        else:
            beta=[beta[0],beta[1]]
    
    if opts.xc==None:
        xc=[im.shape[0]/2.,im.shape[1]/2.]
    else:
        xc=map(float,opts.xc.split(','))
        #corect for position if using only a region of the image
        xc[0]-=extent[0]
        xc[1]-=extent[2]
    if opts.max_pos: xc=shapelets.img.maxPos(im)

    nmax=opts.nmax.split(',')
    if len(nmax)==1:
        nmax=[int(nmax[0])+1,int(nmax[0])+1]
    else:
        nmax=[int(nmax[0])+1,int(nmax[1])+1]
    
    if opts.mode.startswith('pol'):
        beta=(beta[0]+beta[1])/2.
        nmax=nmax[0]
        print "beta0: (%f)\tcentroid: (%f,%f)\tnmax: %i"%(beta,xc[0],xc[1],nmax-1)

        r0,th0=shapelets.shapelet.polarArray(xc,im.shape)

        #plot: data, model, residual: model-data, coeffs
        p.subplot(221)
        p.title('Image')
        p.imshow(im)
        p.colorbar()
        
        p.subplot(222)
        p.title('Model')
        bvals=shapelets.decomp.genPolarBasisMatrix(beta,nmax,r0,th0)
        coeffs=shapelets.decomp.solveCoeffs(bvals,im)
        mdl=n.abs(shapelets.img.constructModel(bvals,coeffs,xc,im.shape))
        p.imshow(mdl)
        p.text(xc[1],xc[0],'+')
        p.colorbar()
        
        p.subplot(223)
        p.title('Residual')
        res=im-mdl
        p.imshow(res)
        p.colorbar()
        
        p.subplot(224)
        p.title('Coefficents')
        cimR=shapelets.img.polarCoeffImg(coeffs.real,nmax)
        cimI=shapelets.img.polarCoeffImg(coeffs.imag,nmax)
        cimI=n.fliplr(cimI)
        cim=n.concatenate((cimR,cimI),axis=1)
        p.pcolor(cim)
        p.colorbar()

        ofn=opts.ofn
        print 'Writing to file:',ofn
        shapelets.fileio.writeLageurreCoeffs(ofn,coeffs,xc,im.shape,beta,nmax,info=ifn)
        
        p.show()
    else:
        print "beta0: (%f,%f)\tcentroid: (%f,%f)\tnmax: (%i,%i)"%(beta[0],beta[1],xc[0],xc[1],nmax[0]-1,nmax[1]-1)

        #plot: data, model, residual: model-data, coeffs
        p.subplot(221)
        p.title('Image')
        p.imshow(im)
        p.colorbar()
        
        p.subplot(222)
        p.title('Model')
        rx=n.array(range(0,im.shape[0]),dtype=float)-xc[0]
        ry=n.array(range(0,im.shape[1]),dtype=float)-xc[1]
        bvals=shapelets.decomp.genBasisMatrix(beta,nmax,rx,ry)
        coeffs=shapelets.decomp.solveCoeffs(bvals,im)
        mdl=shapelets.img.constructModel(bvals,coeffs,xc,im.shape)
        p.imshow(mdl)
        p.text(xc[1],xc[0],'+')
        p.colorbar()
        
        p.subplot(223)
        p.title('Residual')
        res=im-mdl
        p.imshow(res)
        p.colorbar()

        p.subplot(224)
        p.title('Coefficents')
        sqCoeffs=n.reshape(coeffs,nmax)
        p.pcolor(sqCoeffs)
        p.colorbar()
        
        ofn=opts.ofn
        print 'Writing to file:',ofn
        shapelets.fileio.writeHermiteCoeffs(ofn,coeffs,xc,im.shape,beta,nmax,info=ifn)
        
        p.show()

