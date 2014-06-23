"""
Image manipulation functions
"""

import numpy as np
#shapelet functions
import decomp, shapelet

def selPxRange(im,extent):
    """Select out the subimage within the extent region (xmin,xmax,ymin,ymax)"""
    return im[extent[0]:extent[1],extent[2]:extent[3]]

def flux(im):
    """Total flux of the image (M00)"""
    return im.sum()

def centroid(im,region=None):
    """Compute centroid of image (M10/M00,M01/M00) (x,y), region: use this region to compute centroid"""
    if not(region is None):
        im=im[region[0]:region[1],region[2]:region[3]]
        offset=[region[0],region[2]]
    else: offset=[0,0]
    #shift minimum to zero
    im=im-im.min()
    m00 = im.sum()
    m10=np.multiply(np.arange(im.shape[1]),im)
    m01=np.multiply(np.arange(im.shape[0]),np.rot90(im))
    m10=m10.sum()
    m01=m01.sum()
    return [m10/m00+offset[0], m01/m00+offset[1]]

def maxPos(im, region=None):
    """Return the position of the maximum in the image, region: use this region to determine the maximum"""
    if not(region is None):
        im=im[region[0]:region[1],region[2]:region[3]]
        offset=[region[0],region[2]]
    else: offset=[0,0]
    #shift minimum to zero
    im=im-im.min()
    maxpos=np.argwhere(im==np.max(im))[0]
    return [maxpos[0]+offset[0],maxpos[1]+offset[1]]

def estimateNoiseMap(im,region=None,masks=None,sigma=3.,tol=.01,maxiter=None):
    """Generate a noise map based on background pixels by iteratively clipping noise above a set sigma level
    until the variation between iterations is within the tolerance or the maximum number of iterations is reached.
    If region is set then the noise is computed for a region and applied to the entire map.
    Masks can be included to ignore portions of the image."""
    im=np.ma.array(im)
    if region is None:
        if not (masks is None):
            for m in masks:
                im[m[0]:m[1],m[2]:m[3]]=np.ma.masked
        mean0=np.mean(im)
        median0=np.median(im)
        mode0=2.5*median0-1.5*mean0
        std0=np.std(im)
        conv=False
        niter=0
        if maxiter==0:conv=True #compute the noise on the unclipped image
        while not conv:
            print niter
            im=np.ma.masked_greater(im,sigma*np.abs(mode0))
            im=np.ma.masked_less(im,-1*sigma*np.abs(mode0))
            if np.abs(np.std(im)-std0)/std0 < tol: conv=True
            elif np.ma.count_masked(im)>im.size*.5: conv=True
            else:
                std0=np.std(im)
                mode0=2.5*np.median(im)-1.5*np.mean(im)
            niter+=1
            if not(maxiter is None) and niter==maxiter: break
        noisemap=np.ones((im.shape[0],im.shape[1]))*np.std(im)
        return noisemap
    else:
        im_region=selPxRange(im,region)
        std0=np.std(im_region)
        mean0=np.mean(im_region)
        noisemap=np.ones_like(im)
        noisemap=np.random.normal(mean0,std0,(im.shape[0],im.shape[1]))
        return noisemap

def constructModel(bvals,coeffs,xc,size):
    """Construct a model image based on the basis functions values, centroid position xc, and coeffs on an image
    with size dimensions
    """
    model_img=np.dot(bvals,coeffs)
    model_img=np.reshape(model_img,size)
    return model_img

def polarCoeffImg(coeffs,nmax):
    """Return 2D array of coeffs for Laguerre components for plotting
    """
    im=np.zeros((nmax*2,nmax))
    cnt=0
    for nn in range(nmax):
        for mm in np.arange(-1*nn,nn+1):
            if nn%2==0 and mm%2==0:
                im[mm+nmax-1,nn]=coeffs[cnt]
                cnt+=1
            elif nn%2==1 and mm%2==1:
                im[mm+nmax-1,nn]=coeffs[cnt]
                cnt+=1
    return im

def xc2radec(xc,hdr,offset=[0.,0.]):
    """Return the RA,DEC position for a centroid (x,y) pair based on FITS header
    offset: x,y offset from full image
    """
    ra=(hdr['raPix']-(offset[0]+xc[0]))*hdr['dra']+hdr['ra']
    dec=(hdr['decPix']-(offset[1]+xc[1]))*hdr['ddec']+hdr['dec']
    return ra,dec

def beta2size(beta,hdr):
    """Convert a beta pixel size to celestial size
    """
    if type(beta)==list:
        return [np.abs(beta[0]*hdr['dra']),np.abs(beta[1]*hdr['ddec'])]
    else:
        return [np.abs(beta*hdr['dra']),np.abs(beta*hdr['ddec'])]

