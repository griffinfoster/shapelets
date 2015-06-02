"""
Functions for shapelet-based convolution
"""

import sys
import numpy as np
from scipy.misc import factorial

#TODO: polar version
def psfMatrix(g,gamma,alpha,beta,lmax,mmax,nmax,mode='hermite'):
    """Compute the PSF matrix as defined in Refregier and Bacon 2003 Section 3.2
    g: shapelet coefficents of the PSF
    gamma, alpha, beta: scale factors (float)
    nmax: number of coefficients used to represent the convolved image, 2 element list, if 1 element/integer assume it is square
    mmax: number of coefficients used to represent the unconvolved image, 2 element list, if 1 element/integer assume it is square
    lmax: number of coefficients used to represent the PSF, 2 element list, if 1 element/integer assume it is square
    mode: hermite or laguerre
    """
    if mode.startswith('herm'):
        if type(nmax)==int: nmax=[nmax,nmax]
        if type(mmax)==int: nmax=[mmax,mmax]
        if type(lmax)==int: nmax=[lmax,lmax]
        C=generate2dClmnTensor(gamma,alpha,beta,lmax,mmax,nmax) #compute convolution tensor [Refregier and Bacon 2003 eq. 6-11]
        return np.reshape( np.tensordot(C,g,axes=[[0,1],[0,1]]), (C.shape[2]*C.shape[3], C.shape[4]*C.shape[5]) ) #return convolution tensor x g [Refregier and Bacon 2003 section 3.2]

def generateLlmnTensor(a,b,c,lmax,mmax,nmax):
    """Generate the recursive relation L tensor [Refregier and Bacon 2003 eq. 11]
    a,b,c: scale factors (float)
    lmax,mmax,nmax: maximum values for l,m,n (integer)
    """
    
    L=np.zeros((lmax+1,mmax+1,nmax+1))
    #base cases l+m+n=0
    L[0,0,0]=1
    #base cases l+m+n=2
    L[0,1,1]=2.*b*c
    L[1,0,1]=2.*a*c
    L[1,1,0]=2.*a*b
    L[0,0,2]=-2.+2.*(c**2.)
    L[0,2,0]=-2.+2.*(b**2.)
    L[2,0,0]=-2.+2.*(a**2.)

    for n in np.arange(nmax):
        for m in np.arange(mmax):
            for l in np.arange(lmax):
                if l+m+n>2:
                    L[l+1,m,n]=2.*l*((a**2.)-1.)*L[l-1,m,n] + 2.*m*a*b*L[l,m-1,n]           + 2.*n*a*c*L[l,m,n-1]
                    L[l,m+1,n]=2.*l*a*b*L[l-1,m,n]          + 2.*m*((b**2.)-1.)*L[l,m-1,n]  + 2.*n*b*c*L[l,m,n-1]
                    L[l,m,n+1]=2.*l*a*c*L[l-1,m,n]          + 2.*m*b*c*L[l,m-1,n]           + 2.*n*((c**2.)-1.)*L[l,m,n-1]
    return L[:lmax,:mmax,:nmax]

def generateBlmnTensor(a,b,c,lmax,mmax,nmax):
    """Generate the B tensor [Refregier and Bacon 2003 eq. 9]
    a,b,c: scale factors (float)
    lmax,mmax,nmax: maximum values for l,m,n (integer)
    """
    nu=(a**(-2.)+b**(-2.)+c**(-2.))**(-0.5)
    L=generateLlmnTensor(np.sqrt(2.)*nu/a,np.sqrt(2.)*nu/b,np.sqrt(2.)*nu/c,lmax,mmax,nmax)

    lidx, midx, nidx = np.indices((lmax,mmax,nmax))
    B=nu*(((2.**(lidx+midx+nidx-1.))*np.sqrt(np.pi)*factorial(midx)*factorial(nidx)*factorial(lidx)*a*b*c)**(-0.5))*L

    ##slower for loop
    #B=np.zeros((lmax,mmax,nmax))
    #for n in np.arange(nmax):
    #    for m in np.arange(mmax):
    #        for l in np.arange(lmax):
    #            B[l,m,n]=nu*(((2.**(l+m+n-1.))*np.sqrt(np.pi)*factorial(m)*factorial(n)*factorial(l)*a*b*c)**(-0.5))*L[l,m,n]
    return B

def generate1dClmnTensor(a,b,c,lmax,mmax,nmax):
    """Generate the 1-dimensional C tensor [Refregier and Bacon 2003 eq. 7]
    a,b,c: scale factors (float)
    lmax,mmax,nmax: maximum values for l,m,n (integer)
    """
    B=generateBlmnTensor(1./a, 1./b, 1./c, lmax,mmax,nmax)

    lidx, midx, nidx = np.indices((lmax,mmax,nmax))
    C=np.sqrt(2.*np.pi)*((-1.)**nidx)*((1j)**(nidx+midx+lidx))*B

    ##slower for loop
    #C=np.zeros((lmax,mmax,nmax), dtype=np.complex)
    #for n in np.arange(nmax):
    #    for m in np.arange(mmax):
    #        for l in np.arange(lmax):
    #            C[l,m,n]=np.sqrt(2.*np.pi)*((-1.)**n)*((1j)**(n+m+l))*B[l,m,n]
    return C

def generate2dClmnTensor(a,b,c,lmax,mmax,nmax):
    """Generate the 2-dimensional C tensor [Refregier and Bacon 2003 eq. 6]
    a,b,c: scale factors (float)
    lmax,mmax,nmax: maximum values for l,m,n (2 element integer lists)
    """
    C1=generate1dClmnTensor(a,b,c,lmax[0],mmax[0],nmax[0])
    C2=generate1dClmnTensor(a,b,c,lmax[1],mmax[1],nmax[1])
    C=np.kron(C1,C2)
    return np.reshape(C,(lmax[0],lmax[1],mmax[0],mmax[1],nmax[0],nmax[1]))

if __name__ == "__main__":

    print '============================================'
    print 'Testing convolution module:'
    print '============================================'
    tc=0
    te=0

    #generateLlmnTensor():
    tc+=1
    try:
        L=generateLlmnTensor(1.5,2.5,3.5,6,7,8)
        print L.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #generateBlmnTensor():
    tc+=1
    try:
        B=generateBlmnTensor(1.5,2.5,3.5,6,7,8)
        print B.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1
    
    #generate1dClmnTensor():
    tc+=1
    try:
        C1=generate1dClmnTensor(1.5,2.5,3.5,6,7,8)
        print C1.shape
        C2=generate1dClmnTensor(1.5,2.5,3.5,5,4,9)
        print C2.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1
    
    #generate2dClmnTensor():
    tc+=1
    try:
        C=generate2dClmnTensor(1.5,2.5,3.5,[6,5],[7,4],[8,9])
        print C.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #psfMatrix()
    try:
        g=np.random.rand(6,5)
        P=psfMatrix(g,1.5,2.5,3.5,[6,5],[7,4],[8,9],mode='hermite')
        print P.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    Pinv=np.linalg.pinv(P)
    print Pinv.shape

    print '============================================'
    print '%i of %i tests succeeded'%(tc-te,tc)
    print '============================================'

