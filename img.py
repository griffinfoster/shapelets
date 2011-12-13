"""
Image manipulation functions

ToDo:
"""

import numpy as n

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

def estimateNoiseMap(im,masks=None,sigma=3.,tol=.01,maxiter=None):
    """Generate a noise map based on background pixels by iteratively clipping noise above a set sigma level
    until the variation between iterations is with the tolerance or the maximum number of iterations is reached.
    Masks can be included to ignore portions of the image."""
    im=n.ma.array(im)
    if not (masks is None):
        for m in masks:
            im[m[0]:m[1],m[2]:m[3]]=n.ma.masked
    mean0=n.mean(im)
    median0=n.median(im)
    mode0=2.5*median0-1.5*mean0
    std0=n.std(im)
    conv=False
    niter=0
    while not conv:
        im=n.ma.masked_greater(im,sigma*mode0)
        im=n.ma.masked_less(im,-1*sigma*mode0)
        if n.abs(n.std(im)-std0)/std0 < tol: conv=True
        else:
            std0=n.std(im)
            mode0=2.5*n.median(im)-1.5*n.mean(im)
        niter+=1
        if not(maxiter is None) and niter==maxiter: break
    noisemap=n.ones_like(im)
    noisemap=n.random.normal(n.mean(im),n.std(im),(im.shape[0],im.shape[1]))
    return noisemap
    #return n.ma.filled(im,0) #masked array

def constructHermiteModel(n,coeffs,beta,xc,size):
    """Construct a model image based on the beta value, centroid position xc, Hermite basis functions of order n and coeffs on an image
    with size pixels
    """
    
