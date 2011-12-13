"""
Functions for decomposing an image into basis functions using a Least Squares Estimator(LSE)

Implementation of "Statistics in Theory and Practice", Lupton, Ch. 11
"""

import numpy as n

def initBeta(im,frac=.25,nmax=5):
    """Initial starting point for Beta, uses size of image to set limits, initial beta is set to the min beta
    frac: fraction of a pixel to use as the minimum size
    nmax: maximum decomposition order
    beta_max = theta_max / (nmax+1)**.5
    beta_min = theta_min * (nmax+1)**.5
    """
    beta_max=max(im.shape[0]/((nmax+1.)**.5),im.shape[1]/((nmax+1.)**.5))
    return min(frac*((nmax+1.)**.5),beta_max)

def solveCoeffs(m,im):
    """Solve for the basis function coefficents Theta^ using a Least Squares Esimation (Lupton pg. 84)
    theta^ = (m^T * m)^-1 * im
    n: number of pixels in the image
    k: number of basis functions (in all dimensions)
    m: n x k matrix of the basis functions for each pixel
    im: image of size n pixels
    returns theta_hat, a size k array of basis function coefficents"""
    im_flat=n.reshape(im,im.shape[0]*im.shape[1],1)
    mTm=n.dot(m.T,m)                    #matrix multiply m with it's transpose
    mTm_inv=n.linalg.inv(mTm)           #invert the k x k matrix
    mTm_inv_mT=n.dot(mTm_inv,m.T)       #matrix multiply the result with the transpose of m
    theta_hat=n.dot(mTm_inv_mT,im_flat) #compute the coefficents for the basis functions
    return theta_hat

if __name__ == "__main__":
    imgx=50
    imgy=60
    test_im=n.arange(imgx*imgy)
    test_im=n.reshape(test_im,(imgx,imgy))

    nx=5
    ny=6
    test_m=n.arange(nx*ny*imgx*imgy)
    test_m=n.reshape(test_m,(imgx*imgy,nx*ny))
    
    coeffs=solveCoeffs(test_m,test_im)
    print coeffs
    print coeffs.shape
    
