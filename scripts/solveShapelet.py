#!/usr/bin/env python
"""
Solve for shapelet coefficients based on beta, xc, phi, and n_max
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches
import shapelets
import pywcs

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
    o.add_option('--centroid', dest='centroid', action="store_true", default=False,
        help='Use the centroid position instead of max intensity position')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save the figure, requires filename')
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    im0,hdr=shapelets.fileio.readFITS(ifn,hdr=True)
    extent=[0,im0.shape[0],0,im0.shape[1]]

    if not (opts.region is None):
        extent=map(int, opts.region.split(','))
        im=shapelets.img.selPxRange(im0, [extent[2],extent[3],extent[0],extent[1]]) #numpy axis flip)
    else:
        im=im0

    #noise map
    if opts.nregion is None:
        #sample the entire image for noise estimation
        mean,std=shapelets.img.estimateNoise(im0,mode='sample')
        nm=shapelets.img.makeNoiseMap(im.shape,mean,std)
    else:
        #use a specific region for noise estimation
        nextent=map(int, opts.nregion.split(','))
        mean,std=shapelets.img.estimateNoise(shapelets.img.selPxRange(im0,[nextent[2],nextent[3],nextent[0],nextent[1]]),mode='basic')
        nm=shapelets.img.makeNoiseMap(im.shape,mean,std)

    #select initial beta, phi, and xc
    if opts.beta==None:
        beta0,phi0,nmax0=shapelets.decomp.initParams(im,mode='fit',hdr=hdr)
    else:
        beta0=map(float,opts.beta.split(','))
        phi0=float(opts.phi)
        if len(beta0)==1:
            beta0=[beta0[0],beta0[0]]
        else:
            beta0=[beta0[1],beta0[0]] #input to numpy flip

    if opts.centroid:
        xc=shapelets.img.centroid(im)
    elif opts.xc==None:
        xc=shapelets.img.maxPos(im)
    else:
        xc=map(float,opts.xc.split(','))
        xc=[xc[1],xc[0]] #input to numpy flip

    nmax=opts.nmax.split(',')
    if len(nmax)==1:
        nmax=[int(nmax[0])+1,int(nmax[0])+1]
    else:
        nmax=[int(nmax[1])+1,int(nmax[0])+1] #input to numpy flip

    print 'Using beta: (%f,%f) :: \tphi: %f radians :: \tcentre: x,y=(%f,%f) :: \tnmax: (%i,%i)'%(beta0[1],beta0[0],phi0,xc[1],xc[0],nmax[1]-1,nmax[0]-1)

    #determine (RA,dec) coordinates for centroid position
    #TODO: this is correct for when the FITS header is delta RA<0 and delta Dec>0, this may need to be generalized
    if extent is None:
        radec=hdr['wcs'].wcs_pix2sky(np.array([ [xc[1]+1,xc[0]+1] ]),1)[0] #unit: degrees, FITS conventions: first pixel is (1,1)
    else:
        radec=hdr['wcs'].wcs_pix2sky(np.array([ [xc[1]+extent[0]+1,im0.shape[0]-(extent[2]+xc[0])] ]),1)[0] #unit: degrees, FITS conventions: first pixel is (1,1)

    print 'Centroid RA: %f (deg) Dec: %f (deg)'%(radec[0],radec[1])

    if opts.mode.startswith('pol'):
        r0,th0=shapelets.shapelet.polarArray(xc,im.shape)

        #plot: data, model, residual: model-data, coeffs
        fig = plt.figure()
        ax = fig.add_subplot(221)
        plt.title('Image')
        plt.imshow(im)
        e=matplotlib.patches.Ellipse(xy=[xc[1],xc[0]],width=2.*beta0[1],height=2.*beta0[0],angle=(180.*phi0/np.pi)) #numpy to matplotlib flip
        e.set_clip_box(ax.bbox)
        e.set_alpha(0.3)
        e.set_facecolor('black')
        ax.add_artist(e)
        plt.text(xc[1],xc[0],'+',horizontalalignment='center',verticalalignment='center') #numpy to matplotlib flip
        plt.colorbar()
        plt.xlabel('X/RA')
        plt.ylabel('Y/Dec')
        
        plt.subplot(222)
        plt.title('Model')
        bvals=shapelets.decomp.genPolarBasisMatrix(beta0,nmax,phi0,r0,th0)
        coeffs=shapelets.decomp.solveCoeffs(bvals,im)
        mdl=np.abs(shapelets.img.constructModel(bvals,coeffs,im.shape))
        plt.imshow(mdl)
        plt.text(xc[1],xc[0],'+',horizontalalignment='center',verticalalignment='center') #numpy to matplotlib flip
        plt.colorbar()
        plt.xlabel('X/RA')
        plt.ylabel('Y/Dec')
        
        plt.subplot(223)
        plt.title('Residual')
        res=im-mdl
        plt.imshow(res)
        plt.colorbar()
        plt.xlabel('X/RA')
        plt.ylabel('Y/Dec')
        
        plt.subplot(224)
        plt.title('Coefficients')
        cimR=shapelets.img.polarCoeffImg(coeffs.real,nmax)
        cimI=shapelets.img.polarCoeffImg(coeffs.imag,nmax)
        cimI=np.fliplr(cimI)
        cim=np.concatenate((cimR,cimI),axis=1)
        #plt.pcolor(cim)
        plt.imshow(cim,interpolation='nearest',origin='lower')
        plt.colorbar()

        ofn=opts.ofn
        print 'Writing to file:',ofn
        shapelets.fileio.writeLageurreCoeffs(ofn,coeffs,xc,im.shape,beta0,phi0,nmax,info=ifn,pos=[radec[0],radec[1],hdr['dra'],hdr['ddec']])
        
    elif opts.mode.startswith('cart'):

        #plot: data, model, residual: model-data, coeffs
        fig = plt.figure()
        ax = fig.add_subplot(221)
        plt.title('Image')
        plt.imshow(im)
        e=matplotlib.patches.Ellipse(xy=[xc[1],xc[0]],width=2.*beta0[1],height=2.*beta0[0],angle=(180.*phi0/np.pi)) #numpy to matplotlib flip
        e.set_clip_box(ax.bbox)
        e.set_alpha(0.3)
        e.set_facecolor('black')
        ax.add_artist(e)
        plt.text(xc[1],xc[0],'+',horizontalalignment='center',verticalalignment='center') #numpy to matplotlib flip
        plt.colorbar()
        plt.xlabel('X/RA')
        plt.ylabel('Y/Dec')
        
        plt.subplot(222)
        plt.title('Model')
        ry=np.array(range(0,im.shape[0]),dtype=float)-xc[0]
        rx=np.array(range(0,im.shape[1]),dtype=float)-xc[1]
        yy,xx=shapelets.shapelet.xy2Grid(ry,rx)
        bvals=shapelets.decomp.genBasisMatrix(beta0,nmax,phi0,yy,xx)
        coeffs=shapelets.decomp.solveCoeffs(bvals,im)
        mdl=shapelets.img.constructModel(bvals,coeffs,im.shape)
        plt.imshow(mdl)
        plt.text(xc[1],xc[0],'+',horizontalalignment='center',verticalalignment='center') #numpy to matplotlib flip
        plt.colorbar()
        plt.xlabel('X/RA')
        plt.ylabel('Y/Dec')
        
        plt.subplot(223)
        plt.title('Residual')
        res=im-mdl
        plt.imshow(res)
        plt.colorbar()
        plt.xlabel('X/RA')
        plt.ylabel('Y/Dec')

        plt.subplot(224)
        plt.title('Coefficients')
        sqCoeffs=np.reshape(coeffs,nmax)
        #plt.pcolor(sqCoeffs)
        plt.imshow(sqCoeffs,interpolation='nearest',origin='lower')
        plt.colorbar()
        
        ofn=opts.ofn
        print 'Writing to file:',ofn
        shapelets.fileio.writeHermiteCoeffs(ofn,coeffs,xc,im.shape,beta0,phi0,nmax,info=ifn,pos=[radec[0],radec[1],hdr['dra'],hdr['ddec']])
        
    if not (opts.savefig is None):
        plt.savefig(opts.savefig)
    else: plt.show()

