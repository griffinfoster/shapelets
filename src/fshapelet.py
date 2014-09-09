# -*- coding: utf-8 -*-
import numpy as np
from scipy.special import eval_genlaguerre
from scipy.misc import factorial

from numexpr import evaluate as neval

PI = np.pi

def polar_basis_L(r, n0, m0, beta=1.):
    """Polar dimensional basis function based on Laguerre polynomials of characteristic size beta
    phs: additional phase factor, used in the Fourier Transform"""
    m_a = abs(m0)
    n_m = (n0 - m_a)/2
    n_p = (n0 + m_a)/2

    # 2π  is missing here..
    norm = (-1.)**n_m/(beta**(m_a + 1))*(float(factorial(n_m))/factorial(n_p))**.5

    lag = eval_genlaguerre(n_m, m_a, (r/beta)**2)
    return neval('norm*lag*r**m_a*exp(-.5*(r/beta)**2)')

def genPolarBasisMatrix(beta, nmax,r,th):
    """Generate the n x k matrix of basis functions(k) for each pixel(n)
    beta: characteristic size of the shaeplets
    nmax: maximum decomposition order
    r: radius matrix of n pixels
    th: theta matrix of n pixels
    """
    # Precompute the values for angular basis
    # as the complex exp computations are rather costly
    # now Y_{m=-1} can be accessed as Y[-1] or Y[(2*nmax-1) -1]
    Y0 = neval('exp(-1.j*th)')
    N0, N1 = Y0.shape[0], Y0.shape[1]
    Y_vec = []
    for mm in range(nmax) + range(-1*nmax+1, 0):
        Y_vec.append(neval('Y0**mm'))

    bvals=np.empty((N0, N1, pvect_len(nmax)), dtype='complex128')
    L_vec_tmp = [0]*nmax # initializing an empty list
    
    pos = 0
    for nn in range(nmax):
        for mm in range(-1*nn,nn+1):
            if (nn%2==0 and abs(mm)%2==0) or (nn%2==1 and abs(mm)%2==1):
                Ym = wrap_idx(Y_vec, mm)
                # Using the fact that L_{n,-m} = L_{n, m}
                if mm <= 0:
                    Lnm = polar_basis_L(r,nn,mm, beta=beta)
                    L_vec_tmp[abs(mm)] = Lnm
                else:
                    Lnm = L_vec_tmp[mm]
                bvals[:,:,pos] = (neval('Ym*Lnm'))
                pos += 1

    return bvals.reshape((N0*N1,-1))

def wrap_idx(vect, idx):
    # this part reproduces negative indexing
    if idx < 0:
        idx = len(vect) + idx
    return vect[idx]


def pvect_len(nmax):
    return nmax*(nmax+1)/2


def polarArray(xc,size,rot=0.):
    """Return arrays of shape 'size' with radius and theta values centered on xc
    rot: radians in which to rotate the shapelet
    """
    Y, X = np.indices(size, dtype='float64')
    X -= xc[0]
    Y -= xc[1]

    r  = neval('sqrt(X**2 + Y**2)')
    th = neval('arctan2(Y,X)+rot')
    return r, th
