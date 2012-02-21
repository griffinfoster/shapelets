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
    for bid,bf in enumerate(bfs):
        corr+=coeffs[bid]*shapelet.computeBasis2dAtom(bf,uu,vv)
    return n.reshape(corr,dshape)

def computeLaguerreUV(bfs,coeffs,uu,vv):
    """Compute the correlation value for an array of U,V postions based on a set of Laguerre 
    shapelet basis functions and coefficients"""
    dshape=uu.shape
    uu=uu.flatten()
    vv=vv.flatten()
    corr=n.zeros(uu.shape,dtype=complex)
    r0=n.sqrt(n.square(uu) + n.square(vv))
    th0=n.arctan2(vv,uu)
    for bid,bf in enumerate(bfs):
        corr+=coeffs[bid]*shapelet.computeBasisPolarAtom(bf,r0,th0)
    return n.reshape(corr,dshape)

