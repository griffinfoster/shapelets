"""
Functions for Shapelet related operations

ToDo:
    computeBasis: 
        add centroid input for offsetting from the center of an image2D
        add coeff input
    applyCoeff:
        input: coeffs, basis functions
        return basis functions * coeffs
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

def basis2d(n0,n1):
    """2d dimensionless Cartesian basis function"""
    b=hermite2d(n0,n1)
    b[0]*=((2**n0)*(n.pi**(.5))*scipy.factorial(n0))**(-.5)
    exp0=lambda x: b[0](x) * n.exp(-.5*(x**2))
    b[1]*=((2**n1)*(n.pi**(.5))*scipy.factorial(n1))**(-.5)
    exp1=lambda x: b[1](x) * n.exp(-.5*(x**2))
    return [exp0,exp1]

def dimBasis2d(n0,n1,beta=1.):
    """2d dimensional Cartesian basis function of characteristic size beta"""
    b=hermite2d(n0,n1)
    b[0]*=(beta**(-.5))*((2**n0)*(n.pi**(.5))*scipy.factorial(n0))**(-.5)
    exp0=lambda x: b[0](x/beta) * n.exp(-.5*((x/beta)**2))
    b[1]*=(beta**(-.5))*((2**n1)*(n.pi**(.5))*scipy.factorial(n1))**(-.5)
    exp1=lambda x: b[1](x/beta) * n.exp(-.5*((x/beta)**2))
    return [exp0,exp1]

def computeBasis2d(b,rx,ry):
    """Compute the values of a 2D Basis function b in range (rx,ry)"""
    return n.outer(b[0](rx),b[1](ry))

def computeBasis2dAtom(b,x,y):
    """Compute a the basis function in the position (x,y)"""
    return b[0]([x])*b[1]([y])

if __name__ == "__main__":
    xmax=5
    ymax=5
    
    nx=50
    ny=50
    #range is -5/5
    rx=n.array(range(-nx,nx),dtype=float)/10.
    ry=n.array(range(-ny,ny),dtype=float)/10.
    extent=(-5,5,-5,5)

    for y in range(xmax):
        for x in range(ymax):
            b=basis2d(x,y,beta=1.)
            bval=computeBasis2d(b,rx,ry)
            p.subplot(xmax,ymax,(y*ymax+x+1))
            p.imshow(bval,extent=extent,cmap=p.get_cmap('Blues'))
            p.title('(%i,%i)'%(x,y))
    p.show()

