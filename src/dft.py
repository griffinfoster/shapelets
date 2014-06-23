"""
Functions for performing direct Fourier transforms
"""

import numpy as np

def dft2d(d,x,y,k,l,norm=1.):
    """compute the 2d DFT for position (k,l) based on d(x,y)
    norm: normalization factor"""
    return np.sum(d*np.exp(-2.*np.pi*1j*((x*k) + (y*l))))/norm

def idft2d(d,x,y,k,l,norm=1.):
    """compute the 2d Inverse DFT for position (k,l) based on d(x,y)
    norm: normalization factor"""
    return dft2d(d.conjugate(),x,y,k,l,norm=norm).conjugate()

def computeUV(mdl,xx,yy,uu,vv):
    """compute the correlation value for an array of U,V postions based on a model image and 
    the xx,yy coordianates"""
    dshape=uu.shape
    uu=uu.flatten()
    vv=vv.flatten()
    corr=np.zeros(uu.shape,dtype=complex)
    for ind in range(uu.shape[0]):
        corr[ind]=idft2d(mdl,xx,yy,uu[ind],vv[ind],norm=mdl.size)
    return np.reshape(corr,dshape)

