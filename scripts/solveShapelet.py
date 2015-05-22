#!/usr/bin/env python
"""
Solve for shapelet coefficients based on beta, xc, phi, and n_max
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches
import shapelets

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_IMAGE')
    o.set_description(__doc__)
    o.add_option('-r', '--region', dest='region', default=None,
        help='Region of image to decompose into shapelets, (xmin,xmax,ymin,ymax), default: None')
    o.add_option('-N', '--noise_region', dest='nregion', default=None,
        help='Region of image to use to create a noise map, if set to None the entire image is used, this is not used in the script, (xmin,xmax,ymin,ymax), default: None')
    o.add_option('-x', '--xc', dest='xc', default=None,
        help='set a x,y pixel position for shapelet center, if using a region it is based on the relative position, default: centroid of image/region')
    o.add_option('-m', '--mode', dest='mode', default='cart',
        help='Set the shapelet mode, cartesian or polar, default: cartesian')
    o.add_option('-n', '--nmax', dest='nmax', default='5',
        help='Size of coefficient dimensions, can be two values i.e. \'4,5\', default: 5')
    o.add_option('-b', '--beta', dest='beta', default=None,
        help='Beta value, can be two values i.e. \'25.0,30.5\', default: None, guess is made based on Gaussian fit')
    o.add_option('-p','--phi', dest='phi', default=0., type='float',
        help='Rotation angle (radians), only used when beta is manually input, default: 0')
    o.add_option('-o', '--outfile', dest='ofn', default='shapeletCoeffs.pkl',
        help='Coefficients output filename, default: shapeletCoeffs.pkl')
    o.add_option('--max', dest='max_pos', action="store_true", default=False,
        help='Override centroid position to be the position of max intensity')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save the figure, requires filename')
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    im0=shapelets.fileio.readFITS(ifn)
    extent=[0,im0.shape[0],0,im0.shape[1]]
    if not (opts.region is None):
        extent=map(int, opts.region.split(','))
        im=shapelets.img.selPxRange(im0,extent)
    else:
        im=im0

    #noise map
    if opts.nregion is None:
        #use the image region for noise estimation
        mean,std=shapelets.img.estimateNoise(im,mode='basic')
        nm=shapelets.img.makeNoiseMap(im.shape,mean,std)
    else:
        #use a specific region for noise estimation
        nextent=map(int, opts.nregion.split(','))
        mean,std=shapelets.img.estimateNoise(shapelets.img.selPxRange(im0,nextent),mode='basic')
        nm=shapelets.img.makeNoiseMap(im.shape,mean,std)

    #select initial beta, phi, and xc
    if opts.beta==None:
        beta0,phi0=shapelets.decomp.initBetaPhi(im,mode='fit')
    else:
        beta0=map(float,opts.beta.split(','))
        phi0=float(opts.phi)
        if len(beta0)==1:
            beta0=[beta0[0],beta0[0]]
        else:
            beta0=[beta0[0],beta0[1]]
    
    if opts.max_pos:
        xc=shapelets.img.maxPos(im)
    elif opts.xc==None:
        xc=shapelets.img.centroid(im)
    else:
        xc=map(float,opts.xc.split(','))
        #correct for position if using only a region of the image
        #xc[0]-=extent[0]
        #xc[1]-=extent[2]

    nmax=opts.nmax.split(',')
    if len(nmax)==1:
        nmax=[int(nmax[0])+1,int(nmax[0])+1]
    else:
        nmax=[int(nmax[0])+1,int(nmax[1])+1]

    print 'Using beta: (%f,%f) :: \tphi: %f radians :: \tcentre: x,y=(%f,%f) :: \tnmax: (%i,%i)'%(beta0[0],beta0[1],phi0,xc[0],xc[1],nmax[0]-1,nmax[1]-1)

    if opts.mode.startswith('pol'):
        r0,th0=shapelets.shapelet.polarArray(xc,im.shape)

        #plot: data, model, residual: model-data, coeffs
        fig = plt.figure()
        ax = fig.add_subplot(221)
        plt.title('Image')
        plt.imshow(im)
        e=matplotlib.patches.Ellipse(xy=xc,width=2.*np.max(beta0),height=2.*np.min(beta0),angle=(180.*phi0/np.pi))
        e.set_clip_box(ax.bbox)
        e.set_alpha(0.3)
        e.set_facecolor('black')
        ax.add_artist(e)
        plt.text(xc[0],xc[1],'+',horizontalalignment='center',verticalalignment='center')
        plt.colorbar()
        
        plt.subplot(222)
        plt.title('Model')
        bvals=shapelets.decomp.genPolarBasisMatrix(beta0,nmax,phi0,r0,th0)
        coeffs=shapelets.decomp.solveCoeffs(bvals,im)
        mdl=np.abs(shapelets.img.constructModel(bvals,coeffs,im.shape))
        plt.imshow(mdl)
        plt.colorbar()
        
        plt.subplot(223)
        plt.title('Residual')
        res=im-mdl
        plt.imshow(res)
        plt.colorbar()
        
        plt.subplot(224)
        plt.title('Coefficents')
        cimR=shapelets.img.polarCoeffImg(coeffs.real,nmax)
        cimI=shapelets.img.polarCoeffImg(coeffs.imag,nmax)
        cimI=np.fliplr(cimI)
        cim=np.concatenate((cimR,cimI),axis=1)
        plt.pcolor(cim)
        plt.colorbar()

        ofn=opts.ofn
        print 'Writing to file:',ofn
        shapelets.fileio.writeLageurreCoeffs(ofn,coeffs,xc,im.shape,beta0,phi0,nmax,info=ifn)
        
    else:

        #rx=np.array(range(0,im.shape[0]),dtype=float)-xc[0]
        #ry=np.array(range(0,im.shape[1]),dtype=float)-xc[1]
        #xx,yy=shapelets.shapelet.xy2Grid(rx,ry)
        #bvals=shapelets.decomp.genBasisMatrix(beta0,nmax,phi0,xx,yy)
        #fullImg=np.zeros((im.shape[0]*nmax[1],im.shape[1]*nmax[0]))
        #cnt=0
        #for n0 in np.arange(nmax[0]):
        #    for n1 in np.arange(nmax[1]):
        #        fullImg[n0*im.shape[0]:(n0+1)*im.shape[0],n1*im.shape[1]:(n1+1)*im.shape[1]]=np.reshape(bvals[:,cnt],im.shape)
        #        cnt+=1
        #plt.imshow(fullImg)
        #plt.show()

        #plot: data, model, residual: model-data, coeffs
        fig = plt.figure()
        ax = fig.add_subplot(221)
        plt.title('Image')
        plt.imshow(im)
        e=matplotlib.patches.Ellipse(xy=xc,width=2.*np.max(beta0),height=2.*np.min(beta0),angle=(180.*phi0/np.pi))
        e.set_clip_box(ax.bbox)
        e.set_alpha(0.3)
        e.set_facecolor('black')
        ax.add_artist(e)
        plt.text(xc[0],xc[1],'+',horizontalalignment='center',verticalalignment='center')
        plt.colorbar()
        
        plt.subplot(222)
        plt.title('Model')
        rx=np.array(range(0,im.shape[0]),dtype=float)-xc[0]
        ry=np.array(range(0,im.shape[1]),dtype=float)-xc[1]
        xx,yy=shapelets.shapelet.xy2Grid(rx,ry)
        bvals=shapelets.decomp.genBasisMatrix(beta0,nmax,phi0,xx,yy)
        print bvals.shape
        coeffs=shapelets.decomp.solveCoeffs(bvals,im)
        mdl=shapelets.img.constructModel(bvals,coeffs,im.shape)
        plt.imshow(mdl)
        plt.text(xc[0],xc[1],'+',horizontalalignment='center',verticalalignment='center')
        plt.colorbar()
        
        plt.subplot(223)
        plt.title('Residual')
        res=im-mdl
        plt.imshow(res)
        plt.colorbar()

        plt.subplot(224)
        plt.title('Coefficents')
        sqCoeffs=np.reshape(coeffs,nmax)
        plt.pcolor(sqCoeffs)
        plt.colorbar()
        
        ofn=opts.ofn
        print 'Writing to file:',ofn
        shapelets.fileio.writeHermiteCoeffs(ofn,coeffs,xc,im.shape,beta0,nmax,phi0,info=ifn)
        
    if not (opts.savefig is None):
        plt.savefig(opts.savefig)
    else: plt.show()
