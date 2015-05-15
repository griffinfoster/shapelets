"""
Image manipulation functions
"""

import sys
import numpy as np
#shapelet functions
import decomp

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
            print 'iteration:', niter
            im=np.ma.masked_greater(im,sigma*np.abs(mode0))
            im=np.ma.masked_less(im,-1*sigma*np.abs(mode0))
            if np.abs(np.std(im)-std0)/std0 < tol: conv=True
            elif np.ma.count_masked(im)>im.size*.5: conv=True
            else:
                std0=np.std(im)
                mode0=2.5*np.median(im)-1.5*np.mean(im)
            niter+=1
            if not(maxiter is None) and niter==maxiter: break
        #noisemap=np.ones((im.shape[0],im.shape[1]))*np.std(im)
        noisemap=np.random.normal(np.mean(im),np.std(im),(im.shape[0],im.shape[1]))
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

def beta2size(beta,hdr=None,dra=1.,ddec=1.):
    """Convert a beta pixel size to celestial size
    requires either FITS header (hdr) or the delta RA (dra) and delta DEC (ddec)
    """
    if type(beta)==list:
        if hdr is None:
            return [np.abs(beta[0]*dra),np.abs(beta[1]*ddec)]
        else:
            return [np.abs(beta[0]*hdr['dra']),np.abs(beta[1]*hdr['ddec'])]
    else:
        if hdr is None:
            return [np.abs(beta*dra),np.abs(beta*ddec)]
        else:
            return [np.abs(beta*hdr['dra']),np.abs(beta*hdr['ddec'])]

if __name__ == "__main__":

    print '============================================'
    print 'Testing img module:'
    print '============================================'
    tc=0
    te=0

    import fileio
    im,hdr=fileio.readFITS('../data/N6251_test.fits',hdr=True)
    
    #selPxRange(im,extent):
    tc+=1
    try:
        subim=selPxRange(im,[170,335,195,309])
        print subim.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #flux(im):
    tc+=1
    try:
        print flux(im)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #centroid(im,region=None):
    tc+=1
    try:
        print centroid(im,region=[170,335,195,309])
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #maxPos(im, region=None):
    tc+=1
    try:
        print maxPos(im,region=[170,335,195,309])
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #estimateNoiseMap(im,region=None,masks=None,sigma=3.,tol=.01,maxiter=None):
    tc+=1
    try:
        nm=estimateNoiseMap(im)
        print nm.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #constructModel(bvals,coeffs,xc,size):
    tc+=1
    try:
        shapeDict=fileio.readLageurreCoeffs('../data/testHermite.pkl')
        rx=np.array(range(0,shapeDict['size'][0]),dtype=float)-shapeDict['xc'][0]
        ry=np.array(range(0,shapeDict['size'][1]),dtype=float)-shapeDict['xc'][1]
        bvals=decomp.genBasisMatrix(shapeDict['beta'],[shapeDict['norder'],shapeDict['norder']],shapeDict['phi'],rx,ry)
        mdl=constructModel(bvals,shapeDict['coeffs'],shapeDict['xc'],shapeDict['size'])
        print mdl.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #polarCoeffImg(coeffs,nmax):
    tc+=1
    try:
        shapeDict=fileio.readLageurreCoeffs('../data/testLageurre.pkl')
        cim=polarCoeffImg(shapeDict['coeffs'].real,shapeDict['norder'][0])
        print cim
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #xc2radec(xc,hdr,offset=[0.,0.]):
    tc+=1
    try:
        print xc2radec([10.,5.],hdr)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #beta2size
    tc+=1
    try:
        print beta2size([4.5,6.],dra=3.,ddec=4.)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    print '============================================'
    print '%i of %i tests succeeded'%(tc-te,tc)
    print '============================================'

