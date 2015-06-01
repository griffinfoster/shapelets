"""
Image manipulation functions
"""

import sys
import numpy as np
#shapelet functions
import decomp

def selPxRange(im,extent):
    """Select out the subimage within the extent region (xmin,xmax,ymin,ymax)
    Note about x,y def'n: in numpy the top left corner is (0,0) and the first index index increases from top to bottom, the second index increases from left to right, so this is a -90 degree rotation to the normal x-y plane
    """
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
    m00=np.sum(im)
    m01=np.sum(np.sum(im,axis=1)*np.arange(im.shape[0]))
    m10=np.sum(np.sum(im,axis=0)*np.arange(im.shape[1]))
    return [m01/m00+offset[0], m10/m00+offset[1]]

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

def makeNoiseMap(shape,mean=0.,std=1.):
    """Return a noise map with a given shape, mean and std of Gaussian noise
    """
    return np.random.normal(mean,std,shape)

#TODO: auto noise estimation perhaps with edge detection and masking out regions? needs to be worked on, for know we assume a good noise region is known or just sample the entire image
def estimateNoise(im,mode='basic'):
    """
    Estimate the Gaussian noise statistics for an image
    modes:
        basic: compute the mean and std from an input image
        sample: randomly sample a fraction pixels, works fine for sparse images
    """
    if mode.startswith('basic'):
        mean=np.mean(im)
        std=np.std(im)
        print 'Estimated noise: \tmean: %f \tstd: %f'%(mean,std)
        return mean,std
    if mode.startswith('sample'):
        sample=np.random.choice(im.flatten(), size=int(im.size*0.1), replace=True) #randomly sample npix*f pixels
        mean=np.mean(sample)
        std=np.std(sample)
        print 'Estimated noise: \tmean: %f \tstd: %f'%(mean,std)
        return mean,std

def constructModel(bvals,coeffs,size):
    """Construct a model image based on the basis functions values, and coeffs on an image
    with size dimensions
    """
    #return np.reshape(bvals[:,0],size)
    model_img=np.dot(bvals,coeffs)
    model_img=np.reshape(model_img,size)
    return model_img

def polarCoeffImg(coeffs,nmax):
    """Return 2D array of coeffs for Laguerre components for plotting
    """
    im=np.zeros((nmax[0]*2,nmax[0]))
    cnt=0
    for nn in range(nmax[0]):
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
    Notation note: dec is top/bottom direction, ra is left/right direction
    """
    ra=(hdr['raPix']-(offset[1]+xc[1]))*hdr['dra']+hdr['ra']
    dec=(hdr['decPix']-(offset[0]+xc[0]))*hdr['ddec']+hdr['dec']
    return ra,dec

#This has been superseded by measure.beta_pix2angle()
#def beta2size(beta,hdr=None,dra=1.,ddec=1.):
#    """Convert a beta pixel size to celestial size
#    requires either FITS header (hdr) or the delta RA (dra) and delta DEC (ddec)
#    """
#    if type(beta)==list:
#        if hdr is None:
#            return [np.abs(beta[1]*dra),np.abs(beta[0]*ddec)]
#        else:
#            return [np.abs(beta[1]*hdr['dra']),np.abs(beta[0]*hdr['ddec'])]
#    else:
#        if hdr is None:
#            return [np.abs(beta*dra),np.abs(beta*ddec)]
#        else:
#            return [np.abs(beta*hdr['dra']),np.abs(beta*hdr['ddec'])]

if __name__ == "__main__":

    print '============================================'
    print 'Testing img module:'
    print '============================================'
    tc=0
    te=0

    import fileio,shapelet
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

    #estimateNoise(im,mode='basic'):
    tc+=1
    try:
        mean,std=estimateNoise(im,mode='basic')
        nm=makeNoiseMap(im.shape,mean,std)
        print nm.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #constructModel(bvals,coeffs,size):
    tc+=1
    try:
        shapeDict=fileio.readLageurreCoeffs('../data/testHermite.pkl')
        rx=np.array(range(0,shapeDict['size'][0]),dtype=float)-shapeDict['xc'][0]
        ry=np.array(range(0,shapeDict['size'][1]),dtype=float)-shapeDict['xc'][1]
        xx,yy=shapelet.xy2Grid(rx,ry)
        bvals=decomp.genBasisMatrix(shapeDict['beta'],shapeDict['norder'],shapeDict['phi'],xx,yy)
        mdl=constructModel(bvals,shapeDict['coeffs'],shapeDict['size'])
        print mdl.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #polarCoeffImg(coeffs,nmax):
    tc+=1
    try:
        shapeDict=fileio.readLageurreCoeffs('../data/testLageurre.pkl')
        cim=polarCoeffImg(shapeDict['coeffs'].real,shapeDict['norder'])
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

    ##beta2size
    #tc+=1
    #try:
    #    print beta2size([4.5,6.],dra=3.,ddec=4.)
    #except:
    #    print 'Test failed (%i):'%tc, sys.exc_info()[0]
    #    te+=1

    print '============================================'
    print '%i of %i tests succeeded'%(tc-te,tc)
    print '============================================'

