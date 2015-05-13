"""
Functions for decomposing an image into basis functions using a Least Squares Estimator(LSE)

Implementation of "Statistics in Theory and Practice", Lupton, Ch. 11
"""

import numpy as np
import shapelet,img
import sys

#TODO; check these initBeta functions
def initBeta(im,frac=.25,nmax=5):
    """Initial starting point for Beta, uses size of image to set limits, initial beta is set to the min beta
    frac: fraction of a pixel to use as the minimum size
    nmax: maximum decomposition order
    beta_max = theta_max / (nmax+1)**.5
    beta_min = theta_min * (nmax+1)**.5
    """
    beta_max=max(im.shape[0]/((nmax+1.)**.5),im.shape[1]/((nmax+1.)**.5))
    #print im.shape[0]/((nmax+1.)**.5),im.shape[1]/((nmax+1.)**.5)
    beta0=min(frac*((nmax+1.)**.5),beta_max)
    return [beta0,beta0]

def initBeta2(im,frac=.2):
    """Initial starting point for Beta, uses size of image
    frac: fraction of image to use as the initial beta
    """
    return [frac*im.shape[0],frac*im.shape[1]]

def genPolarBasisMatrix(beta,nmax,phi,r,th):
    """Generate the n x k matrix of basis functions(k) for each pixel(n)
    beta: characteristic size of the shapelets (b_major, b_minor)
    nmax: maximum decomposition order
    phi: rotation angle
    r: radius matrix of n pixels
    th: theta matrix of n pixels
    """
    bvals=[]
    for nn in range(nmax):
        for mm in np.arange(-1*nn,nn+1):
            if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                bf=shapelet.polarDimBasis(nn,mm,beta=beta,phi=phi)
                bvals.append(shapelet.computeBasisPolar(bf,r,th).flatten())
    bm=np.array(bvals)
    return bm.transpose()

def genBasisMatrix(beta,nmax,phi,rx,ry):
    """Generate the n x k matrix of basis functions(k) for each pixel(n)
    nmax: maximum decompisition order
    beta: characteristic size of the shapelet
    phi: rotation angle
    rx: range of x values to evaluate basis functions
    ry: range of y values to evaluate basis functions
    """
    bvals=[]
    for x in range(nmax[0]):
        for y in range(nmax[1]):
            bf=shapelet.dimBasis2d(x,y,beta=beta,phi=phi)
            bvals.append(shapelet.computeBasis2d(bf,rx,ry).flatten())
    bm=np.array(bvals)
    return bm.transpose()

def solveCoeffs(m,im):
    """Solve for the basis function coefficients Theta^ using a Least Squares Esimation (Lupton pg. 84)
    theta^ = (m^T * m)^-1 * im
    n: number of pixels in the image
    k: number of basis functions (in all dimensions)
    m: n x k matrix of the basis functions for each pixel
    im: image of size n pixels
    returns theta_hat, a size k array of basis function coefficents"""
    #im_flat=np.reshape(im,im.shape[0]*im.shape[1],1)
    im_flat=im.flatten()
    mTm=np.dot(m.T,m)                    #matrix multiply m with it's transpose
    mTm_inv=np.linalg.inv(mTm)           #invert the k x k matrix
    mTm_inv_mT=np.dot(mTm_inv,m.T)       #matrix multiply the result with the transpose of m
    theta_hat=np.dot(mTm_inv_mT,im_flat) #compute the coefficents for the basis functions
    return theta_hat

def chi2PolarFunc(params,nmax,im,nm):
    """Function which is to be minimized in the chi^2 analysis for Polar shapelets
    params = [beta0, beta1, phi, xc, yc]
        beta0: characteristic size of shapelets, fit parameter
        beta1: characteristic size of shapelets, fit parameter
        phi: rotation angle of shapelets, fit parameter
        xc: x centroid of shapelets, fit parameter
        yc: y centroid of shapelets, fit parameter
    nmax: number of coefficents to use in the Laguerre polynomials
    im: observed image
    nm: noise map
    """
    beta0=params[0]
    beta1=params[1]
    phi=params[2]
    xc=params[3]
    yc=params[4]
    if beta0<0.:
        print 'warning: beta going negative, setting to 0.0'
        beta0=0.
    if beta1<0.:
        print 'warning: beta going negative, setting to 0.0'
        beta1=0.
    print 'beta: (%f,%f)\t phi: %f\txc: (%f,%f)'%(beta0,beta1,phi,xc,yc)

    size=im.shape
    r,th=shapelet.polarArray([xc,yc],size)
    bvals=genPolarBasisMatrix([beta0,beta1],nmax,phi,r,th)
    coeffs=solveCoeffs(bvals,im)
    mdl=np.abs(img.constructModel(bvals,coeffs,[xc,yc],size))
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2betaPolarFunc(params,xc,yc,r,th,nmax,im,nm):
    """Function which is to be minimized in the chi^2 analysis for Polar shapelets
    params = [beta0, beta1, phi]
        beta0: characteristic size of shapelets, fit parameter
        beta1: characteristic size of shapelets, fit parameter
        phi: rotation angle of shapelets, fit parameter
    xc: x centroid of shapelets
    yc: y centroid of shapelets
    r: radius from centroid, array of im.shape
    th: angle from centroid, array of im.shape
    nmax: number of coefficents to use in the Laguerre polynomials
    im: observed image
    nm: noise map
    """
    beta0=params[0]
    beta1=params[1]
    phi=params[2]
    if beta0<0.:
        print 'warning: beta going negative, setting to 0.0'
        beta0=0.
    if beta1<0.:
        print 'warning: beta going negative, setting to 0.0'
        beta1=0.
    print 'beta: (%f,%f)\t phi: %f'%(beta0,beta1,phi)

    size=im.shape
    bvals=genPolarBasisMatrix([beta0,beta1],nmax,phi,r,th)
    coeffs=solveCoeffs(bvals,im)
    mdl=np.abs(img.constructModel(bvals,coeffs,[xc,yc],size))
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2nmaxPolarFunc(params,im,nm,beta0,beta1,phi,xc):
    """
    params = [nmax]
        nmax: number of coefficents
    im: observed image
    nm: noise map
    beta0: fit beta value
    beta1: fit beta value
    phi: rotation angle
    xc: fit centroid position
    """
    #print params
    nmax=params[0]
    size=im.shape
    r,th=shapelet.polarArray(xc,size)
    bvals=genPolarBasisMatrix([beta0,beta1],nmax,phi,r,th)
    coeffs=solveCoeffs(bvals,im)
    mdl=np.abs(img.constructModel(bvals,coeffs,xc,size))
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2Func(params,nmax,im,nm):
    """Function which is to be minimized in the chi^2 analysis
    params = [betaX, betaY, phi, xc, yc]
        betaX: characteristic size of shapelets, fit parameter
        betaY: characteristic size of shapelets, fit parameter
        phi: rotation angle of shapelets, fit parameter
        xc: x centroid of shapelets, fit parameter
        yc: y centroid of shapelets, fit parameter
    nmax: number of coefficents to use in x,y
    im: observed image
    nm: noise map
    """
    betaX=params[0]
    betaY=params[1]
    phi=params[2]
    xc=params[3]
    yc=params[4]
    if betaX<0.:
        print 'warning: betaX going negative, setting to 0.0'
        betaX=0.
    if betaY<0.:
        print 'warning: betaY going negative, setting to 0.0'
        betaY=0.
    print 'beta: (%f,%f)\tphi: %f\txc: (%f,%f)'%(betaX,betaY,phi,xc,yc)
    
    size=im.shape
    #shift the (0,0) point to the centroid
    rx=np.array(range(0,size[0]),dtype=float)-xc
    ry=np.array(range(0,size[1]),dtype=float)-yc

    bvals=genBasisMatrix([betaX,betaY],nmax,phi,rx,ry)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,[xc,yc],size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2betaFunc(params,xc,yc,nmax,im,nm):
    """Function which is to be minimized in the chi^2 analysis
    params = [betaX, betaY, phi]
        betaX: characteristic size of shapelets, fit parameter
        betaY: characteristic size of shapelets, fit parameter
        phi: rotation angle of shapelets, fit parameter
    xc: x centroid of shapelets
    yc: y centroid of shapelets
    nmax: number of coefficents to use in x,y
    im: observed image
    nm: noise map
    """
    betaX=params[0]
    betaY=params[1]
    phi=params[2]
    if betaX<0.:
        print 'warning: betaX going negative, setting to 0.0'
        betaX=0.
    if betaY<0.:
        print 'warning: betaY going negative, setting to 0.0'
        betaY=0.
    print 'beta: (%f,%f)\tphi: %f'%(betaX,betaY,phi)
    size=im.shape
    #shift the (0,0) point to the centroid
    rx=np.array(range(0,size[0]),dtype=float)-xc
    ry=np.array(range(0,size[1]),dtype=float)-yc

    bvals=genBasisMatrix([betaX,betaY],nmax,phi,rx,ry)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,[xc,yc],size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2xcFunc(params,beta0,beta1,phi,nmax,im,nm):
    """Function which is to be minimized in the chi^2 analysis
    params = [xc, yc]
        xc: x centroid of shapelets, fit parameter
        yc: y centroid of shapelets, fit parameter
    beta0: characteristic size of shapelet
    beta1: characteristic size of shapelet
    phi: rotation angle of shapelets
    nmax: number of coefficents to use in x,y
    im: observed image
    nm: noise map
    """
    xc=params[0]
    yc=params[1]
    print xc,yc
    size=im.shape
    #shift the (0,0) point to the centroid
    rx=np.array(range(0,size[0]),dtype=float)-xc
    ry=np.array(range(0,size[1]),dtype=float)-yc

    bvals=genBasisMatrix([beta0,beta1],nmax,phi,rx,ry)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,[xc,yc],size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2nmaxFunc(params,im,nm,beta0,beta1,phi,xc):
    """
    params = [nmaxX,nmaxY]
        nmax: number of coefficents to use in x,y
    im: observed image
    nm: noise map
    beta0: characteristic size of shapelet
    beta1: characteristic size of shapelet
    phi: rotation angle of shapelets
    xc: fit centroid position
    """
    print params
    nmax=params
    size=im.shape
    #shift the (0,0) point to the centroid
    rx=np.array(range(0,size[0]),dtype=float)-xc[0]
    ry=np.array(range(0,size[1]),dtype=float)-xc[1]

    bvals=genBasisMatrix([beta0,beta1],nmax,phi,rx,ry)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,xc,size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

if __name__ == "__main__":

    print '============================================'
    print 'Testing decomp module:'
    print '============================================'
    import fileio
    tc=0
    te=0
    
    im,hdr=fileio.readFITS('../data/N6251_test.fits',hdr=True)
    subim=img.selPxRange(im,[170,335,195,309])
    
    #guess beta
    #initBeta(im,frac=.25,nmax=5):
    #initBeta2(im,frac=.2):
    tc+=1
    try:
        print initBeta(subim,frac=.25,nmax=5)
        print initBeta2(subim,frac=.2)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #genPolarBasisMatrix(beta,nmax,r,th):
    #solveCoeffs(m,im):
    tc+=1
    try:
        beta0=initBeta2(subim,frac=.2)
        xc=img.maxPos(subim)
        r0,th0=shapelet.polarArray(xc,subim.shape)
        mPolar=genPolarBasisMatrix(beta0,5,0.,r0,th0)
        coeffs=solveCoeffs(mPolar,subim)
        print coeffs
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #genBasisMatrix(beta,nmax,rx,ry):
    #solveCoeffs(m,im):
    tc+=1
    try:
        beta0=initBeta2(subim,frac=.2)
        xc=img.maxPos(subim)
        rx=np.array(range(0,subim.shape[0]),dtype=float)-xc[0]
        ry=np.array(range(0,subim.shape[1]),dtype=float)-xc[1]
        mCart=genBasisMatrix(beta0,[5,5],0.,rx,ry)
        coeffs=solveCoeffs(mCart,subim)
        print coeffs
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    beta0=initBeta2(subim,frac=.2)
    xc=img.maxPos(subim)
    nm=img.estimateNoiseMap(im,region=[170,335,195,309])
    nm=img.selPxRange(nm,[170,335,195,309])

    #chi2PolarFunc(params,nmax,im,nm):
    tc+=1
    try:
        func0=chi2PolarFunc([beta0[0],beta0[1],0.,xc[0],xc[1]],5,subim,nm)
        print func0
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #chi2betaPolarFunc(params,xc,yc,r,th,nmax,im,nm):
    tc+=1
    try:
        r0,th0=shapelet.polarArray(xc,subim.shape)
        func0=chi2betaPolarFunc([beta0[0],beta0[1],0.],xc[0],xc[1],r0,th0,5,subim,nm)
        print func0
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #chi2nmaxPolarFunc(params,im,nm,beta,xc):
    tc+=1
    try:
        print chi2nmaxPolarFunc([5],subim,nm,beta0[0],beta0[1],0.,xc)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #chi2Func(params,nmax,im,nm):
    tc+=1
    try:
        func0=chi2Func([beta0[0],beta0[1],0.,xc[0],xc[1]],[5,5],subim,nm)
        print func0
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #chi2betaFunc(params,xc,yc,nmax,im,nm):
    tc+=1
    try:
        func0=chi2betaFunc([beta0[0],beta0[1],0.],xc[0],xc[1],[5,5],subim,nm)
        print func0
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #chi2xcFunc(params,beta0,beta1,phi,nmax,im,nm):
    tc+=1
    try:
        print chi2xcFunc([xc[0],xc[1]],beta0[0],beta0[1],0.,[5,5],subim,nm)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #chi2nmaxFunc(params,im,nm,beta,xc):
    tc+=1
    try:
        print chi2nmaxFunc([5,5],subim,nm,beta0[0],beta0[1],0.,xc)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    print '============================================'
    print '%i of %i tests succeeded'%(tc-te,tc)
    print '============================================'

