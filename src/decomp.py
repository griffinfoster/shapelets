"""
Functions for decomposing an image into basis functions using a Least Squares Estimator(LSE)

Implementation of "Statistics in Theory and Practice", Lupton, Ch. 11
"""

import numpy as np
import shapelet,img

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

def genPolarBasisMatrix(beta,nmax,r,th):
    """Generate the n x k matrix of basis functions(k) for each pixel(n)
    beta: characteristic size of the shaeplets
    nmax: maximum decomposition order
    r: radius matrix of n pixels
    th: theta matrix of n pixels
    """
    bvals=[]
    for nn in range(nmax):
        for mm in np.arange(-1*nn,nn+1):
            if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                bf=shapelet.polarDimBasis(nn,mm,beta=beta)
                bvals.append(shapelet.computeBasisPolar(bf,r,th).flatten())
    bm=np.array(bvals)
    return bm.transpose()

def genBasisMatrix(beta,nmax,rx,ry):
    """Generate the n x k matrix of basis functions(k) for each pixel(n)
    nmax: maximum decompisition order
    beta: characteristic size of the shapelet
    rx: range of x values to evaluate basis functions
    ry: range of y values to evaluate basis functions
    """
    bvals=[]
    for x in range(nmax[0]):
        for y in range(nmax[1]):
            bf=shapelet.dimBasis2d(x,y,beta=beta)
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
    params = [beta, xc, yc]
        beta: characteristic size of shapelets, fit parameter
        xc: x centroid of shapelets, fit parameter
        yc: y centroid of shapelets, fit parameter
    nmax: number of coefficents to use in the Laguerre polynomials
    im: observed image
    nm: noise map
    """
    beta=params[0]
    xc=params[1]
    yc=params[2]
    if beta<0.:
        print 'warning: beta going negative, setting to 0.0'
        beta=0.
    print 'beta: %f\txc: (%f,%f)'%(beta,xc,yc)

    size=im.shape
    r,th=shapelet.polarArray([xc,yc],size)
    bvals=genPolarBasisMatrix(beta,nmax,r,th)
    coeffs=solveCoeffs(bvals,im)
    mdl=np.abs(img.constructModel(bvals,coeffs,[xc,yc],size))
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2betaPolarFunc(params,xc,yc,r,th,nmax,im,nm):
    """Function which is to be minimized in the chi^2 analysis for Polar shapelets
    params = [beta]
        beta: characteristic size of shapelets, fit parameter
    xc: x centroid of shapelets
    yc: y centroid of shapelets
    r: radius from centroid, array of im.shape
    th: angle from centroid, array of im.shape
    nmax: number of coefficents to use in the Laguerre polynomials
    im: observed image
    nm: noise map
    """
    beta=params[0]
    if beta<0.:
        print 'warning: beta going negative, setting to 0.0'
        beta=0.
    print 'beta: %f'%(beta)

    size=im.shape
    bvals=genPolarBasisMatrix(beta,nmax,r,th)
    coeffs=solveCoeffs(bvals,im)
    mdl=np.abs(img.constructModel(bvals,coeffs,[xc,yc],size))
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2nmaxPolarFunc(params,im,nm,beta,xc):
    """
    params = [nmax]
        nmax: number of coefficents
    im: observed image
    nm: noise map
    beta: fit beta value
    xc: fit centroid position
    """
    print params
    nmax=params
    size=im.shape
    r,th=shapelet.polarArray(xc,size)
    bvals=genPolarBasisMatrix(beta,nmax,r,th)
    coeffs=solveCoeffs(bvals,im)
    mdl=np.abs(img.constructModel(bvals,coeffs,xc,size))
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2Func(params,nmax,im,nm):
    """Function which is to be minimized in the chi^2 analysis
    params = [beta, xc, yc]
        beta: characteristic size of shapelets, fit parameter
        xc: x centroid of shapelets, fit parameter
        yc: y centroid of shapelets, fit parameter
    nmax: number of coefficents to use in x,y
    im: observed image
    nm: noise map
    """
    betaX=params[0]
    betaY=params[1]
    xc=params[2]
    yc=params[3]
    if betaX<0.:
        print 'warning: betaX going negative, setting to 0.0'
        betaX=0.
    if betaY<0.:
        print 'warning: betaY going negative, setting to 0.0'
        betaY=0.
    print 'beta: (%f,%f)\txc: (%f,%f)'%(betaX,betaY,xc,yc)
    
    size=im.shape
    #shift the (0,0) point to the centroid
    rx=np.array(range(0,size[0]),dtype=float)-xc
    ry=np.array(range(0,size[1]),dtype=float)-yc

    bvals=genBasisMatrix([betaX,betaY],nmax,rx,ry)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,[xc,yc],size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2betaFunc(params,xc,yc,nmax,im,nm):
    """Function which is to be minimized in the chi^2 analysis
    params = [beta]
        beta: characteristic size of shapelets, fit parameter
    xc: x centroid of shapelets
    yc: y centroid of shapelets
    nmax: number of coefficents to use in x,y
    im: observed image
    nm: noise map
    """
    betaX=params[0]
    betaY=params[1]
    if betaX<0.:
        print 'warning: betaX going negative, setting to 0.0'
        betaX=0.
    if betaY<0.:
        print 'warning: betaY going negative, setting to 0.0'
        betaY=0.
    print 'beta: (%f,%f)'%(betaX,betaY)
    size=im.shape
    #shift the (0,0) point to the centroid
    rx=np.array(range(0,size[0]),dtype=float)-xc
    ry=np.array(range(0,size[1]),dtype=float)-yc

    bvals=genBasisMatrix([betaX,betaY],nmax,rx,ry)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,[xc,yc],size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2xcFunc(params,beta,nmax,im,nm):
    """Function which is to be minimized in the chi^2 analysis
    params = [xc, yc]
        xc: x centroid of shapelets, fit parameter
        yc: y centroid of shapelets, fit parameter
    beta: characteristic size of shapelet
    nmax: number of coefficents to use in x,y
    im: observed image
    nm: noise map
    """
    xc=params[2]
    yc=params[3]
    print xc,yc
    size=im.shape
    #shift the (0,0) point to the centroid
    rx=np.array(range(0,size[0]),dtype=float)-xc
    ry=np.array(range(0,size[1]),dtype=float)-yc

    bvals=genBasisMatrix(beta,nmax,rx,ry)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,[xc,yc],size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2nmaxFunc(params,im,nm,beta,xc):
    """
    params = [nmax]
        nmax: number of coefficents to use in x,y
    im: observed image
    nm: noise map
    beta: fit beta values
    xc: fit centroid position
    """
    print params
    nmax=[params,params]
    size=im.shape
    #shift the (0,0) point to the centroid
    rx=np.array(range(0,size[0]),dtype=float)-xc[0]
    ry=np.array(range(0,size[1]),dtype=float)-xc[1]

    bvals=genBasisMatrix(beta,nmax,rx,ry)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,xc,size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

if __name__ == "__main__":
    print 'testing decomp'
    #call fitting scripts
