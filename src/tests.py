#!/usr/bin/python
# -*- coding: utf-8 -*-
import nose
import numpy as np

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in, assert_less_equal, \
    assert_greater_equal
from numpy.testing import assert_allclose

import fshapelet as fsh
import shapelet as sh
import cshapelet as csh
import decomp


def setup():
    global r, th, xc, beta, PB_DIM
    xc = (100, 130)
    beta = 100.
    PB_DIM = (300, 200)
    r, th = sh.polarArray(xc, PB_DIM)

def test_polar_vectors_numexpr():
    """ Check that the the version with numexpr generates the same output """
    for nn in range(5):
        for mm in range(-1*nn,nn+1):
            if (nn%2==0 and abs(mm)%2==0) or (nn%2==1 and abs(mm)%2==1):
                res0 = fsh.polar_basis_L(r,nn,mm, beta)*np.exp(-1j*th*mm)
                res1 = sh.polarDimBasis(nn,mm, beta)(r,th)
                yield assert_allclose, res0, res1, 1e-7, 1e-9

def test_polar_basis_numexpr():
    nmax = 10
    res0 = fsh.genPolarBasisMatrix(beta, nmax, r, th)
    res1 = decomp.genPolarBasisMatrix(beta, nmax, r, th)
    yield assert_allclose, res0, res1, 1e-7, 1e-9

def test_polar_basis_cython():
    nmax = 10
    res0 = csh.genPolarBasisMatrix(beta, nmax, r, th)
    res1 = decomp.genPolarBasisMatrix(beta, nmax, r, th)
    yield assert_allclose, res0, res1, 1e-7, 1e-9

def test_polarArray():
    xc = (100., 200.)
    res0 = fsh.polarArray(xc, PB_DIM)
    res1 = sh.polarArray(xc, PB_DIM)
    yield assert_allclose, res0, res1, 1e-7, 1e-9





