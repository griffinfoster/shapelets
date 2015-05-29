"""
Measurement operations
"""

#TODO: shear estimator
#TODO: convert beta to angle on the sky, currently they are in pixels
#TODO: R**2 returns negative values
import numpy as np
import scipy.special

def polarIDs(nmax):
    """Given an nmax, convert to (n.m) polar ID pairs
    """
    ids=[]
    for nn in range(nmax[0]):
        for mm in np.arange(-1*nn,nn+1):
            if nn%2==0 and mm%2==0:
                ids.append([nn,mm])
            elif nn%2==1 and mm%2==1:
                ids.append([nn,mm])
    return ids

def cartIDs(nmax):
    """Given an nmax, convert to (n,m) Cartesian ID pairs
    """
    ids=[]
    if type(nmax) is int: nmax=[nmax,nmax]
    for ny in range(nmax[0]):
        for nx in range(nmax[1]):
            ids.append([ny,nx])
    return ids

def flux(coeffs,beta,nmax,mode):
    """Total flux from Cartesian or polar shapelet coefficients [manual eq. 1.14 and 1.23]
    coeffs: shapelet coefficients
    beta: two element beta
    nmax: coefficient limit
    mode: Hermite or Lageurre
    """
    if mode.startswith('hermite'):
        c0=np.reshape(coeffs,nmax)
        flux=0
        for ny in range(nmax[0]):
            for nx in range(nmax[1]):
                if ny%2==0 and nx%2==0:
                    flux+=(2.**(.5*(2-ny-nx)))*np.sqrt(scipy.special.binom(ny,ny/2))*np.sqrt(scipy.special.binom(nx,nx/2))*c0[ny,nx]
        return np.sqrt(np.pi*beta[0]*beta[1])*flux
    elif mode.startswith('lageurre'):
        ids=polarIDs(nmax)
        flux=0
        for cnt,i in enumerate(ids):
            if i[0]%2==0 and i[1]==0: flux+=coeffs[cnt]
        return np.sqrt(4.*np.pi*beta[0]*beta[1])*flux.real
    else:
        print 'Error: Unknown mode'
        return np.nan

def centroid(coeffs,beta,nmax,mode):
    """Centroid position from Cartesian or polar shapelet coefficients [manual eq. 1.15 and 1.24]
    coeffs: shapelet coefficients
    beta: two element beta
    nmax: coefficient limit
    mode: Hermite or Lageurre
    """
    f=flux(coeffs,beta,nmax,mode)
    if mode.startswith('hermite'):
        c0=np.reshape(coeffs,nmax)
        xc=np.array([0.,0.])
        for ny in range(nmax[0]):
            for nx in range(nmax[1]):
                if ny%2==1 and nx%2==0:
                    xc[0]+=np.sqrt(ny+1)*(2**(.5*(2-ny-nx)))*np.sqrt(scipy.special.binom(ny+1,(ny+1)/2))*np.sqrt(scipy.special.binom(nx,nx/2))*c0[ny,nx]
                elif ny%2==0 and nx%2==1:
                    xc[1]+=np.sqrt(nx+1)*(2**(.5*(2-ny-nx)))*np.sqrt(scipy.special.binom(ny,ny/2))*np.sqrt(scipy.special.binom(nx+1,(nx+1)/2))*c0[ny,nx]
        return (1./f)*np.sqrt(np.pi)*beta[0]*beta[1]*xc
    elif mode.startswith('lageurre'):
        ids=polarIDs(nmax)
        xc=0
        for cnt,i in enumerate(ids):
            if i[0]%2==1 and i[1]==1: xc+=np.sqrt(i[0]+1)*coeffs[cnt]
        xc=((np.sqrt(8*np.pi)*beta[0]*beta[1])/f)*xc
        return np.array([xc.real,xc.imag])
    else:
        print 'Error: Unknown mode'
        return np.nan

def quadrupoles(coeffs,beta,nmax,mode='hermite'):
    """Quadrupole moments (J11, J12, J21, J22) from Cartesian shapelet coefficients [manual eq. 1.16 and 1.17]
    coeffs: shapelet coefficients
    beta: two element beta
    nmax: coefficient limit
    mode: only Hermite supported
    """
    jj=np.zeros([2,2])
    f=flux(coeffs,beta,nmax,mode)
    if mode.startswith('hermite'):
        c0=np.reshape(coeffs,nmax)
        xc=np.array([0.,0.])
        for ny in range(nmax[0]):
            for nx in range(nmax[1]):
                if ny%2==1 and nx%2==1:
                    jj[0,1]=np.sqrt(ny+1)*np.sqrt(nx+1)*(2.**(.5*(2.-ny-nx)))*np.sqrt(scipy.special.binom(ny+1,(ny+1)/2))*np.sqrt(scipy.special.binom(nx+1,(nx+1)/2))*c0[ny,nx]
                    jj[1,0]=np.sqrt(ny+1)*np.sqrt(nx+1)*(2.**(.5*(2.-ny-nx)))*np.sqrt(scipy.special.binom(ny+1,(ny+1)/2))*np.sqrt(scipy.special.binom(nx+1,(nx+1)/2))*c0[ny,nx]
                elif ny%2==0 and nx%2==0:
                    jj[0,0]=(2.*ny+1)*(2.**(.5*(2.-ny-nx)))*np.sqrt(scipy.special.binom(ny,ny/2))*np.sqrt(scipy.special.binom(nx,nx/2))*c0[ny,nx]
                    jj[1,1]=(2.*nx+1)*(2.**(.5*(2.-ny-nx)))*np.sqrt(scipy.special.binom(ny,ny/2))*np.sqrt(scipy.special.binom(nx,nx/2))*c0[ny,nx]
    else:
        print 'Error: Unknown mode'
        return np.nan
    return np.sqrt(np.pi*((beta[0]*beta[1])**3.))*(1/f)*jj

def r2size(coeffs,beta,nmax,mode):
    """Object size from Cartesian or polar shapelet coefficients [manual eq. 1.18 and 1.25]
    coeffs: shapelet coefficients
    beta: two element beta
    nmax: coefficient limit
    mode: Hermite or Lageurre
    """
    f=flux(coeffs,beta,nmax,mode)
    if mode.startswith('hermite'):
        c0=np.reshape(coeffs,nmax)
        r2size=0
        for ny in range(nmax[0]):
            for nx in range(nmax[1]):
                if ny%2==0 and nx%2==0:
                    r2size+=(2.**(.5*(4.-ny-nx)))*(1.+ny+nx)*np.sqrt(scipy.special.binom(ny,ny/2))*np.sqrt(scipy.special.binom(nx,nx/2))*c0[ny,nx]
        return np.sqrt(np.pi*((beta[0]*beta[1])**3.))*(1./f)*r2size
    elif mode.startswith('lageurre'):
        ids=polarIDs(nmax)
        r2size=0
        for cnt,i in enumerate(ids):
            if i[0]%2==0 and i[1]==0: r2size+=(i[0]+1)*coeffs[cnt]
        return np.sqrt(16.*np.pi*((beta[0]*beta[1])**3.))*(1./f)*r2size.real
    else:
        print 'Error: Unknown mode'
        return np.nan

def ellipticity(coeffs,beta,nmax,mode='lageurre'):
    """Object ellipticity from polar shapelet coefficients [manual eq. 1.26]
    coeffs: shapelet coefficients
    beta: two element beta
    nmax: coefficient limit
    mode: only Laguerre supported
    """
    f=flux(coeffs,beta,nmax,mode)
    r2=r2size(coeffs,beta,nmax,mode)
    if mode.startswith('lageurre'):
        ids=polarIDs(nmax)
        ee=0.
        for cnt,i in enumerate(ids):
            if i[0]%2==0 and i[1]==2: ee+=np.sqrt(i[0]*(i[0]+1.))*coeffs[cnt]
        return (np.sqrt(16.*np.pi*((beta[0]*beta[1])**3.))/(f*r2))*ee
    else:
        print 'Error: Unknown mode'
        return np.nan

def all(coeffs,beta,nmax,mode):
    """Comppute all measurements for a set of coefficients
    coeffs: shapelet coefficients
    beta: two element beta
    nmax: coefficient limit
    mode: Hermite or Lageurre
    """
    if mode.startswith('hermite'):
        return {'flux': flux(coeffs,beta,nmax,mode),
                'centroid': centroid(coeffs,beta,nmax,mode),
                'quadrupoles': quadrupoles(coeffs,beta,nmax,mode),
                'r2': r2size(coeffs,beta,nmax,mode) }
    elif mode.startswith('lageurre'):
        return {'flux': flux(coeffs,beta,nmax,mode),
                'centroid': centroid(coeffs,beta,nmax,mode),
                'r2': r2size(coeffs,beta,nmax,mode),
                'ellipticity': ellipticity(coeffs,beta,nmax,mode)}
    else:
        print 'Error: Unknown mode'
        return np.nan

#TODO: cartesian to polar coeff transform
#def cart2polar(coeffs,nmax):
#    """Convert coefficients from Hermite Cartesian to Lageurre Polar
#    coeffs: shapelet coefficients
#    nmax: coefficient limit
#    """
#    ids=polarIDs(nmax)
#    pcoeffs=np.zeros(len(ids))
#    for cnt,i in enumerate(ids):
        

if __name__ == "__main__":

    print '============================================'
    print 'Testing measure module:'
    print '============================================'
    import fileio
    import sys
    tc=0
    te=0
    
    #load precomputed shapelet coeffs (polar and cartesian)
    hermDict=fileio.readHermiteCoeffs('../data/testHermite.pkl')
    laDict=fileio.readLageurreCoeffs('../data/testLageurre.pkl')

    #polarIDs(nmax):
    tc+=1
    try:
        print polarIDs(5)
        print cartIDs(5)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #compute flux (cartesian)
    #compute flux (polar)
    tc+=1
    try:
        print flux(hermDict['coeffs'],hermDict['beta'],hermDict['norder'],mode=hermDict['mode'])
        print flux(laDict['coeffs'],laDict['beta'],laDict['norder'][0],mode=laDict['mode'])
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #compute centroid (polar)
    #compute centroid (cartesian)
    tc+=1
    try:
        print centroid(hermDict['coeffs'],hermDict['beta'],hermDict['norder'],mode=hermDict['mode'])
        print centroid(laDict['coeffs'],laDict['beta'],laDict['norder'][0],mode=laDict['mode'])
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #compute quadrupoles (cartesian J11, J12, J21, J22)
    tc+=1
    try:
        print quadrupoles(hermDict['coeffs'],hermDict['beta'],hermDict['norder'])
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1
    
    #compute size (cartesian)
    #compute size (polar)
    tc+=1
    try:
        print r2size(hermDict['coeffs'],hermDict['beta'],hermDict['norder'],mode=hermDict['mode'])
        print r2size(laDict['coeffs'],laDict['beta'],laDict['norder'][0],mode=laDict['mode'])
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #compute ellipticity (polar)
    tc+=1
    try:
        print ellipticity(laDict['coeffs'],laDict['beta'],laDict['norder'][0])
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    tc+=1
    try:
        print all(hermDict['coeffs'],hermDict['beta'],hermDict['norder'],mode=hermDict['mode'])
        print all(laDict['coeffs'],laDict['beta'],laDict['norder'][0],mode=laDict['mode'])
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #print cart2polar(hermDict['coeffs'],hermDict['norder'])

    print '============================================'
    print '%i of %i tests succeeded'%(tc-te,tc)
    print '============================================'

