"""
Functions for Shapelet related operations
"""

import numpy as n
from scipy.misc import factorial
import scipy.special as special

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
    b[0]*=((2**n0)*(n.pi**(.5))*factorial(n0))**(-.5)
    exp0=lambda x: beta[0] * b[0](x) * n.exp(-.5*(x**2))
    b[1]*=((2**n1)*(n.pi**(.5))*factorial(n1))**(-.5)
    exp1=lambda x: beta[1] * b[1](x) * n.exp(-.5*(x**2))
    return [exp0,exp1]

def dimBasis2d(n0,n1,beta=[1.,1.],phs=[1.,1.]):
    """2d dimensional Cartesian basis function of characteristic size beta
    phs: additional phase factor, used in the Fourier Transform"""
    b=hermite2d(n0,n1)
    b[0]*=(beta[0]**(-.5))*(((2**n0)*(n.pi**(.5))*factorial(n0))**(-.5))
    exp0=lambda x: b[0](x/beta[0]) * n.exp(-.5*((x/beta[0])**2)) * phs[0]
    b[1]*=(beta[1]**(-.5))*(((2**n1)*(n.pi**(.5))*factorial(n1))**(-.5))
    exp1=lambda x: b[1](x/beta[1]) * n.exp(-.5*((x/beta[1])**2)) * phs[1]
    return [exp0,exp1]

def polarDimBasis(n0,m0,beta=1.,phs=1.):
    """Polar dimensional basis function based on Laguerre polynomials of characteristic size beta
    phs: additional phase factor, used in the Fourier Transform"""
    b0=laguerre(n0,m0)
    norm=(((-1.)**((n0-n.abs(m0))/2))/(beta**(n.abs(m0)+1)))*((float(factorial(int((n0-n.abs(m0))/2)))/float(factorial(int((n0+n.abs(m0))/2))))**.5)
    exp0=lambda r,th: norm * r**(n.abs(m0)) * b0((r**2.)/(beta**2.)) * n.exp(-.5*(r**2.)/(beta**2.)) * n.exp(-1j*m0*th)
    return exp0

def polarArray(xc,size,rot=0.):
    """Return arrays of shape 'size' with radius and theta values centered on xc
    rot: radians in which to rotate the shapelet
    """
    rx=n.array(range(0,size[1]),dtype=float)-xc[0]
    ry=n.array(range(0,size[0]),dtype=float)-xc[1]
    rx=n.reshape(n.tile(rx,size[0]),(size[0],size[1]))
    ry=n.reshape(n.tile(ry,size[1]),(size[1],size[0]))
    rExp = lambda x,y: n.sqrt(n.square(x) + n.square(y))
    thExp = lambda x,y: n.arctan2(y,x)+rot
    return rExp(rx,ry.T), thExp(rx,ry.T)

def xy2rth(rx,ry):
    """Convert a range of x and y to r,th arrays of shape (len(x),len(y))"""
    rx0=n.reshape(n.tile(rx,len(ry)),(len(ry),len(rx)))
    ry0=n.reshape(n.tile(ry,len(rx)),(len(rx),len(ry)))
    rExp = lambda x,y: n.sqrt(n.square(x) + n.square(y))
    thExp = lambda x,y: n.arctan2(y,x)
    return rExp(rx0,ry0.T), thExp(rx0,ry0.T)

def rth2xy(r,th):
    """Convert r,theta array pair to an x,y pair"""
    x=r*n.cos(th)
    y=r*n.sin(th)
    return x,y

def xyRotate(rx,ry,rot=0.):
    """Apply a rotation(radians) to an set of X,Y coordinates"""
    r0,th0=xy2rth(rx,ry)
    th0+=rot
    return rth2xy(r0,th0)

def computeBasisPolar(b,r,th):
    """Compute the values of a Polar Basis function b over the R and Theta range"""
    return b(r,th)

def computeBasisPolarAtom(b,r,th):
    """Compute the polar basis function b in the position (rad,theta)"""
    return b(r,th)

def computeBasis2d(b,rx,ry):
    """Compute the values of a 2D Basis function b in range (rx,ry)"""
    return n.outer(b[0](rx),b[1](ry))

def computeBasis2dAtom(b,x,y):
    """Compute the basis function b in the position (x,y), x and y can be arrays"""
    return (b[0]([x])*b[1]([y]))[0]

def ftHermiteBasis(beta,nmax):
    """generate a set of Fourier Transformed Hermite basis functions
    nmax: maximum decompisition order
    beta: characteristic size of the shapelet
    """
    bfs=[]
    for x in range(nmax[0]):
        for y in range(nmax[1]):
            bfs.append(dimBasis2d(x,y,beta=[1./beta[0],1./beta[1]],phs=[1j**(x),1j**(y)]))
            #bfs.append(dimBasis2d(x,y,beta=[1./beta[0],1./beta[1]],phs=[1j**(x+1),1j**(y+1)]))
            #bfs.append(dimBasis2d(x,y,beta=[1./beta[0],1./beta[1]]))
    return bfs

def ftLaguerreBasis(beta,nmax):
    """generate a set of Fourier Transformed Laguerre basis functions
    nmax: maximum decompisition order
    beta: characteristic size of the shapelet
    """
    bfs=[]
    for nn in range(nmax):
        for mm in n.arange(-1*nn,nn+1):
            if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                bfs.append(polarDimBasis(nn,mm,beta=(1./beta),phs=(1j**nn * 1j**mm)))
    return bfs

if __name__ == "__main__":
    print 'shapelets main'

