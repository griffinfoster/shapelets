# -*- coding: utf-8 -*-
#cython: profile=True
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

import numpy as np
cimport numpy as np
from scipy.special import eval_genlaguerre
from scipy.misc import factorial

cimport cython
from cython.parallel import parallel, prange
from libc.math cimport exp, cos, sin, abs
import numexpr as ne

from fshapelet import polar_basis_L

PI = np.pi

cdef extern from "complex.h":
    double complex cexp(double complex z) nogil

null_char_array = np.array([[0]], dtype='uint8')

cpdef polar_basis_Y0(double [:,::1] th, double phs=1.,
                            char [:,::1] mask=null_char_array):
    """
    Basically complex exp computations taking into account a mask
    """
    cdef int k, l
    cdef int N0,N1
    cdef char use_mask= mask.size > 1
    N0, N1 = th.shape[0], th.shape[1]

    cdef double complex [:,::1] vect=np.empty((N0, N1), dtype='complex128')*np.nan
    for k in range(N0):
        for l in range(N1):
            if use_mask:
                if mask[k,l]==0:
                    continue
            vect[k,l] = cexp(-1j*th[k,l])

    return vect.base


cpdef genPolarBasisMatrix(double beta, int nmax,r,th, mask=None):
    """Generate the n x k matrix of basis functions(k) for each pixel(n)
    beta: characteristic size of the shaeplets
    nmax: maximum decomposition order
    r: radius matrix of n pixels
    th: theta matrix of n pixels
    """
    cdef int nn, mm, N0, N1, pos
    if mask is None:
        mask = null_char_array
    # Precompute the values for angular basis
    # as the complex exp computations are rather costly
    # now Y_{m=-1} can be accessed as Y[-1] or Y[(2*nmax-1) -1]
    Y0 = polar_basis_Y0(th, mask=mask)
    N0, N1 = Y0.shape[0], Y0.shape[1]
    Y_vec = []
    for mm in range(nmax) + range(-1*nmax+1, 0):
        Y_vec.append(Y0**mm)

    bvals=np.empty((N0*N1, pvect_len(nmax)), dtype='complex128')
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
                bvals[:,pos] = (Ym*Lnm).flatten()
                pos += 1

    return bvals

cdef wrap_idx(vect, int idx):
    # this part reproduces negative indexing
    if idx < 0:
        idx = len(vect) + idx
    return vect[idx]


cdef int pvect_len(int nmax):
    return nmax*(nmax+1)/2
