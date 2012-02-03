"""
Functions for Shapelet related operations
"""

import numpy as n
import pylab as p
import scipy
import scipy.special as special
import math

def hermite2d(n0,n1):
    """Return a n0 x n1 order 2D Hermite polynomial"""
    h0=special.hermite(n0)
    h1=special.hermite(n1)
    return [h0,h1]

def laguerre(n0,m0):
    """Return a generalized Laguerre polynomial L^(|m|)_((n-|m|)/2)(x)"""
    l0=special.genlaguerre(n=(n0-n.abs(m0))/2,alpha=n.abs(m0))
    return l0

def basis2d(n0,n1,beta=[1.,1.]):
    """2d dimensionless Cartesian basis function"""
    b=hermite2d(n0,n1)
    b[0]*=((2**n0)*(n.pi**(.5))*scipy.factorial(n0))**(-.5)
    exp0=lambda x: beta[0] * b[0](x) * n.exp(-.5*(x**2))
    b[1]*=((2**n1)*(n.pi**(.5))*scipy.factorial(n1))**(-.5)
    exp1=lambda x: beta[1] * b[1](x) * n.exp(-.5*(x**2))
    return [exp0,exp1]

def dimBasis2d(n0,n1,beta=[1.,1.]):
    """2d dimensional Cartesian basis function of characteristic size beta"""
    b=hermite2d(n0,n1)
    b[0]*=(beta[0]**(-.5))*(((2**n0)*(n.pi**(.5))*scipy.factorial(n0))**(-.5))
    exp0=lambda x: b[0](x/beta[0]) * n.exp(-.5*((x/beta[0])**2))
    b[1]*=(beta[1]**(-.5))*(((2**n1)*(n.pi**(.5))*scipy.factorial(n1))**(-.5))
    exp1=lambda x: b[1](x/beta[1]) * n.exp(-.5*((x/beta[1])**2))
    return [exp0,exp1]

def polarDimBasis(n0,m0,beta=1.):
    """Polar dimensional basis function based on Laguerre polynomials of characteristic size beta"""
    b0=laguerre(n0,m0)
    norm=(((-1.)**((n0-n.abs(m0))/2))/(beta**(n.abs(m0)+1)))*((float(scipy.factorial(int((n0-n.abs(m0))/2)))/float(scipy.factorial(int((n0+n.abs(m0))/2))))**.5)
    exp0=lambda r,th: norm * r**(n.abs(m0)) * b0((r**2.)/(beta**2.)) * n.exp(-.5*(r**2.)/(beta**2.)) * n.exp(-1j**m0*th)
    return exp0

def computeBasis2d(b,rx,ry):
    """Compute the values of a 2D Basis function b in range (rx,ry)"""
    return n.outer(b[0](rx),b[1](ry))

def computeBasis2dAtom(b,x,y):
    """Compute a the basis function in the position (x,y)"""
    return b[0]([x])*b[1]([y])

if __name__ == "__main__":
    print 'shapelets main'

