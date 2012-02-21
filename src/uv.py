"""
Functions relate to computing UV data from shapelets
"""

import numpy as n
import shapelet

def computeHermiteUV(bfs,coeffs,uu,vv):
    """Compute the correlation value for an array of U,V postions based on a set of Hermite 
    shapelet basis functions and coefficients"""
    dshape=uu.shape
    uu=uu.flatten()
    vv=vv.flatten()
    corr=n.zeros(uu.shape,dtype=complex)
    for ind in range(uu.shape[0]):
        for bid,bf in enumerate(bfs):
            corr[ind]+=coeffs[bid]*shapelet.computeBasis2dAtom(bf,uu[ind],vv[ind])
    return n.reshape(corr,dshape)

def computeLaguerreUV(bfs,coeffs,uu,vv):
    """Compute the correlation value for an array of U,V postions based on a set of Laguerre 
    shapelet basis functions and coefficients"""
    dshape=uu.shape
    uu=uu.flatten()
    vv=vv.flatten()
    corr=n.zeros(uu.shape,dtype=complex)
    for ind in range(uu.shape[0]):
        r0=n.sqrt(n.square(uu[ind]) + n.square(vv[ind]))
        th0=n.arctan2(vv[ind],uu[ind])
        for bid,bf in enumerate(bfs):
            corr[ind]+=coeffs[bid]*shapelet.computeBasisPolarAtom(bf,r0,th0)
    return n.reshape(corr,dshape)

