"""
Functions relate to computing UV data from shapelets
"""

import numpy as np
import shapelet

def computeHermiteUV(bfs,coeffs,uu,vv):
    """Compute the correlation value for an array of U,V postions based on a set of Hermite 
    shapelet basis functions and coefficients
    """
    dshape=uu.shape
    uu=uu.flatten()
    vv=vv.flatten()
    vis=np.zeros(uu.shape,dtype=complex)
    for bid,bf in enumerate(bfs):
        vis+=coeffs[bid]*shapelet.computeBasis2dAtom(bf,uu,vv)
    return np.reshape(vis,dshape)

def computeLaguerreUV(bfs,coeffs,uu,vv):
    """Compute the correlation value for an array of U,V postions based on a set of Laguerre 
    shapelet basis functions and coefficients"""
    dshape=uu.shape
    uu=uu.flatten()
    vv=vv.flatten()
    vis=np.zeros(uu.shape,dtype=complex)
    r0=np.sqrt(np.square(uu) + np.square(vv))
    th0=np.arctan2(vv,uu)
    for bid,bf in enumerate(bfs):
        vis+=coeffs[bid]*shapelet.computeBasisPolarAtom(bf,r0,th0)
    return np.reshape(vis,dshape)

