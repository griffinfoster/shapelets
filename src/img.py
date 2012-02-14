"""
Image manipulation functions
"""

import numpy as n
import pylab as p
import time
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
    m10=n.multiply(n.arange(im.shape[1]),im)
    m01=n.multiply(n.arange(im.shape[0]),n.rot90(im))
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
    maxpos=n.argwhere(im==n.max(im))[0]
    return [maxpos[0]+offset[0],maxpos[1]+offset[1]]

def estimateNoiseMap(im,region=None,masks=None,sigma=3.,tol=.01,maxiter=None):
    """Generate a noise map based on background pixels by iteratively clipping noise above a set sigma level
    until the variation between iterations is within the tolerance or the maximum number of iterations is reached.
    If region is set then the noise is computed for a region and applied to the entire map.
    Masks can be included to ignore portions of the image."""
    im=n.ma.array(im)
    if region is None:
        if not (masks is None):
            for m in masks:
                im[m[0]:m[1],m[2]:m[3]]=n.ma.masked
        mean0=n.mean(im)
        median0=n.median(im)
        mode0=2.5*median0-1.5*mean0
        std0=n.std(im)
        conv=False
        niter=0
        if maxiter==0:conv=True #compute the noise on the unclipped image
        #imHist=[]
        while not conv:
            im=n.ma.masked_greater(im,sigma*n.abs(mode0))
            im=n.ma.masked_less(im,-1*sigma*n.abs(mode0))
            if n.abs(n.std(im)-std0)/std0 < tol: conv=True
            else:
                std0=n.std(im)
                mode0=2.5*n.median(im)-1.5*n.mean(im)
            niter+=1
            #imHist.append(im)
            if not(maxiter is None) and niter==maxiter: break
        noisemap=n.ones_like(im)
        noisemap=n.random.normal(n.mean(im),n.std(im),(im.shape[0],im.shape[1]))
        #for h in range(len(imHist)):
        #    p.subplot(1,len(imHist),h+1)
        #    p.imshow(imHist[h])
        #p.show()
        return noisemap
    else:
        im_region=selPxRange(im,region)
        std0=n.std(im_region)
        mean0=n.mean(im_region)
        noisemap=n.ones_like(im)
        noisemap=n.random.normal(mean0,std0,(im.shape[0],im.shape[1]))
        return noisemap

def constructModel(bvals,coeffs,xc,size):
    """Construct a model image based on the basis functions values, centroid position xc, and coeffs on an image
    with size dimensions
    """
    model_img=n.dot(bvals,coeffs)
    model_img=n.reshape(model_img,size)
    return model_img

def polarCoeffImg(coeffs,nmax):
    """Return 2D array of coeffs for Laguerre components for plotting
    """
    im=n.zeros((nmax*2,nmax))
    cnt=0
    for nn in range(nmax):
        for mm in n.arange(-1*nn,nn+1):
            if nn%2==0 and mm%2==0:
                im[mm+nmax-1,nn]=coeffs[cnt]
                cnt+=1
            elif nn%2==1 and mm%2==1:
                im[mm+nmax-1,nn]=coeffs[cnt]
                cnt+=1
    return im

