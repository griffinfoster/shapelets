"""
Functions for Shapelet related operations
"""

import sys

import numpy as np
from scipy.misc import factorial
from scipy import special

#TODO: hermite 2d, round gaussian?
#TODO: Fourier transform
#########################################################
#def polar2cart():
#    """Convert a set of polar coefficients to Cartesian coefficients [manual eq. 1.27]
#    """

def hermite2d(n0,n1):
    """Return a n0 x n1 order 2D Hermite polynomial"""
    h0=special.hermite(n0)
    h1=special.hermite(n1)
    return [h0,h1]

def laguerre(n0,m0):
    """Return a generalized Laguerre polynomial L^(|m|)_((n-|m|)/2)(x)"""
    l0=special.genlaguerre(n=(n0-np.abs(m0))/2,alpha=np.abs(m0))
    return l0

def rotMatrix(phi):
    """2D Cartesian rotation matrix (radians)"""
    return np.matrix([[np.cos(phi),-1.*np.sin(phi)],[np.sin(phi),np.cos(phi)]])

def basis2d(n0,n1,beta=[1.,1.],phi=0.,fourier=False):
    """2d dimensionless Cartesian basis function
    phi: rotation angle
    fourier: return the Fourier transformed version of the function
    """
    b=hermite2d(n0,n1)
    m=rotMatrix(phi)
    phs=1.
    if fourier:
        beta=[1./beta[0],1./beta[1]]
        phs=1j**(n0+n1)
    b[0]*=((2**n0)*(np.pi**(.5))*factorial(n0))**(-.5)*phs
    exp0=lambda x: beta[0] * b[0](x) * np.exp(-.5*(x**2))
    b[1]*=((2**n1)*(np.pi**(.5))*factorial(n1))**(-.5)*phs
    exp1=lambda x: beta[1] * b[1](x) * np.exp(-.5*(x**2))

    return lambda x,y: exp0(m[0,0]*x+m[0,1]*y)*exp1(m[1,0]*x+m[1,1]*y)

def dimBasis2d(n0,n1,beta=[1.,1.],phi=0.,fourier=False):
    """2d dimensional Cartesian basis function of characteristic size beta
    phi: rotation angle
    fourier: return the Fourier transformed version of the function
    """
    b=hermite2d(n0,n1)
    m=rotMatrix(phi)
    phs=1.
    if fourier:
        beta=[1./beta[0],1./beta[1]]
        phs=1j**(n0+n1)
    b[0]*=(beta[0]**(-.5))*(((2**n0)*(np.pi**(.5))*factorial(n0))**(-.5))*phs
    exp0=lambda x: b[0](x/beta[0]) * np.exp(-.5*((x/beta[0])**2))
    b[1]*=(beta[1]**(-.5))*(((2**n1)*(np.pi**(.5))*factorial(n1))**(-.5))*phs
    exp1=lambda x: b[1](x/beta[1]) * np.exp(-.5*((x/beta[1])**2))

    return lambda x,y: exp0(m[0,0]*x+m[0,1]*y)*exp1(m[1,0]*x+m[1,1]*y)

#TODO: make into an elliptical form?
def polarDimBasis(n0,m0,beta=1.,phi=0.,fourier=False):
    """Polar dimensional basis function based on Laguerre polynomials of characteristic size beta
    phi: rotation angle
    fourier: return the Fourier transformed version of the function
    """
    if len(beta)==1: beta=[beta,beta]
    phs=1.
    if fourier:
        beta=[1./beta[0],1./beta[1]]
        phs=1j**(n0+m0)
    b0=laguerre(n0,m0)
    norm=(((-1.)**((n0-np.abs(m0))/2))/np.sqrt(beta[0]**(np.abs(m0)+1)*beta[1]**(np.abs(m0)+1)))*((float(factorial(int((n0-np.abs(m0))/2)))/float(factorial(int((n0+np.abs(m0))/2))))**.5)*phs
    exp0=lambda r,th: norm * r**(np.abs(m0)) * b0((r**2.)/(beta[0]*beta[1])) * np.exp(-.5*(r**2.)/(beta[0]*beta[1])) * np.exp(-1j*m0*(th+phi))
    return exp0

def polarArray(xc,size,rot=0.):
    """Return arrays of shape 'size' with radius and theta values centered on xc
    rot: radians in which to rotate the shapelet
    """
    rx=np.array(range(0,size[1]),dtype=float)-xc[0]
    ry=np.array(range(0,size[0]),dtype=float)-xc[1]
    rx=np.reshape(np.tile(rx,size[0]),(size[0],size[1]))
    ry=np.reshape(np.tile(ry,size[1]),(size[1],size[0]))
    rExp = lambda x,y: np.sqrt(np.square(x) + np.square(y))
    thExp = lambda x,y: np.arctan2(y,x)+rot
    return rExp(rx,ry.T), thExp(rx,ry.T)

def cartArray(xc,size):
    """Return arrays of shape 'size' with x,y values centered on xc
    """
    rx=np.array(range(0,size[1]),dtype=float)-xc[0]
    ry=np.array(range(0,size[0]),dtype=float)-xc[1]
    rx=np.reshape(np.tile(rx,size[0]),(size[0],size[1]))
    ry=np.reshape(np.tile(ry,size[1]),(size[1],size[0]))
    return rx,ry.T

def xy2Grid(rx,ry):
    """Convert a range of x and y to a grid of shape (len(x),len(y))"""
    rx0=np.reshape(np.tile(rx,len(ry)),(len(ry),len(rx)))
    ry0=np.reshape(np.tile(ry,len(rx)),(len(rx),len(ry)))
    return rx0,ry0.T

def xy2rthGrid(rx,ry):
    """Convert a range of x and y to r,th arrays of shape (len(x),len(y))"""
    rx0=np.reshape(np.tile(rx,len(ry)),(len(ry),len(rx)))
    ry0=np.reshape(np.tile(ry,len(rx)),(len(rx),len(ry)))
    rExp = lambda x,y: np.sqrt(np.square(x) + np.square(y))
    thExp = lambda x,y: np.arctan2(y,x)
    return rExp(rx0,ry0.T), thExp(rx0,ry0.T)

def rth2xy(r,th):
    """Convert r,theta array pair to an x,y pair"""
    x=r*np.cos(th)
    y=r*np.sin(th)
    return x,y

def xyRotate(rx,ry,rot=0.):
    """Apply a rotation(radians) to an set of X,Y coordinates"""
    r0,th0=xy2rthGrid(rx,ry)
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
    return b(rx,ry)

def computeBasis2dAtom(b,x,y):
    """Compute the basis function b in the position (x,y), x and y can be arrays"""
    return b(rx,ry)

if __name__ == "__main__":

    print '============================================'
    print 'Testing shapelets module:'
    print '============================================'
    tc=0
    te=0

    #hermite2d(n0,n1):
    tc+=1
    try:
        h0,h1=hermite2d(3,4)
        print 'hermite:', type(h0), type(h1)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #laguerre(n0,m0):
    tc+=1
    try:
        l0=laguerre(3,3)
        print 'laguerre:', type(l0)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #rotMatrix(phi):
    tc+=1
    try:
        print rotMatrix(np.pi/4.)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #xyRotate(rx,ry,rot=0.):
    tc+=1
    try:
        xp,yp=xyRotate(np.array([1.]),np.array([0.]),rot=np.pi/4.)
        print xp,yp
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #basis2d(n0,n1,beta=[1.,1.],phi=0.):
    tc+=1
    try:
        b=basis2d(3,4,beta=[1.,1.],phi=np.pi/4.)
        print b(2.,3.5)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1
        
    #dimBasis2d(n0,n1,beta=[1.,1.],phi=0.):
    tc+=1
    try:
        b=dimBasis2d(3,4,beta=[1.,1.],phi=np.pi/4.)
        print b(2.,3.5)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1
    
    #polarDimBasis(n0,m0,beta=[1.,1.],phi=0.):
    tc+=1
    try:
        b=polarDimBasis(3,3,beta=[1.,1.],phi=np.pi/4.)
        print b(3.5,np.pi/8.)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1
    
    #polarArray(xc,size,rot=0.):
    tc+=1
    try:
        r,th=polarArray([0.,0.],[15,20],rot=0.)
        print r.shape, th.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #cartArray(xc,size):
    tc+=1
    try:
        x,y=cartArray([6.,7.],[15,20])
        print x.shape, y.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #xy2rthGrid(rx,ry):
    tc+=1
    try:
        r,th=xy2rthGrid(np.arange(10),np.arange(10))
        print r.shape, th.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #rth2xy(r,th):
    tc+=1
    try:
        x,y=rth2xy(np.random.randn(10),2.*np.pi*np.random.rand(10))
        print x.shape,y.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #computeBasisPolar(b,r,th):
    tc+=1
    try:
        r,th=polarArray([0.,0.],[15,20],rot=0.)
        b=polarDimBasis(3,3,beta=[1.,1.],phi=np.pi/4.)
        bval=computeBasisPolar(b,r,th)
        print bval.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #computeBasisPolarAtom(b,r,th):
    tc+=1
    try:
        b=polarDimBasis(3,3,beta=[1.,1.],phi=np.pi/4.)
        bval=computeBasisPolar(b,5.,np.pi/8.)
        print bval
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #computeBasis2d(b,rx,ry):
    tc+=1
    try:
        rx,ry=cartArray([6.,7.],[15,20])
        b=dimBasis2d(3,4,beta=[1.,1.],phi=np.pi/4.)
        bval=computeBasis2d(b,rx,ry)
        print bval.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #computeBasis2dAtom(b,x,y):
    tc+=1
    try:
        b=dimBasis2d(3,4,beta=[1.,1.],phi=np.pi/4.)
        bval=computeBasis2dAtom(b,np.random.randn(10),np.random.randn(10))
        print bval.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    print '============================================'
    print '%i of %i tests succeeded'%(tc-te,tc)
    print '============================================'
