"""
Functions for decomposing an image into basis functions using a Least Squares Estimator(LSE)

Implementation of "Statistics in Theory and Practice", Lupton, Ch. 11
"""

import numpy as np
import shapelet,img
from scipy import optimize

#TODO; a better initBeta function
#def initBeta(im,frac=.25,nmax=5):
#    """Initial starting point for Beta, uses size of image to set limits, initial beta is set to the min beta
#    frac: fraction of a pixel to use as the minimum size
#    nmax: maximum decomposition order
#    beta_max = theta_max / (nmax+1)**.5
#    beta_min = theta_min * (nmax+1)**.5
#    """
#    beta_max=max(im.shape[0]/((nmax+1.)**.5),im.shape[1]/((nmax+1.)**.5))
#    #print im.shape[0]/((nmax+1.)**.5),im.shape[1]/((nmax+1.)**.5)
#    beta0=min(frac*((nmax+1.)**.5),beta_max)
#    return [beta0,beta0]
#
#def initBeta(im,frac=.2):
#    """Initial starting point for Beta, uses size of image
#    frac: fraction of image to use as the initial beta
#    """
#    return [frac*im.shape[0],frac*im.shape[1]]

def ellipticalGaussian2D(x0=0.,y0=0.,sigmax=1.,sigmay=1.,phi=0.,amp=1.,offset=0.):
    """A generalized 2D ellipitical Gaussian function with centre (x0,y0), width (sigmax,sigmay), amplitude amp, offset, and rotation angle phi
    """
    aa=((np.cos(phi))**2.)/(2.*(sigmax**2.)) + ((np.sin(phi))**2.)/(2.*(sigmay**2.))
    bb=(-1.*np.sin(2.*phi))/(4.*(sigmax**2.)) + (np.sin(2.*phi))/(4.*(sigmay**2.))
    cc=((np.sin(phi))**2.)/(2.*(sigmax**2.)) + ((np.cos(phi))**2.)/(2.*(sigmay**2.))
    return lambda x,y: amp * np.exp( -1.*( aa*((x-x0)**2.) + 2.*bb*(x-x0)*(y-y0) + cc*((y-y0)**2.) ) ) + offset

def initGaussian(im):
    """Returns (x, y, sigmax, sigmay, phi, amp, offset)
    the gaussian parameters of a 2D distribution by calculating its
    moments
    borrowed from: http://wiki.scipy.org/Cookbook/FittingData"""
    total = im.sum()
    X, Y = np.indices(im.shape)
    x = (X*im).sum()/total
    y = (Y*im).sum()/total
    col = im[:, int(y)]
    sigmax = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/np.abs(col).sum())
    row = im[int(x), :]
    sigmay = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/np.abs(row).sum())
    offset=np.median(im)
    amp = im.max()-offset
    phi=0.
    return x, y, sigmax, sigmay, phi, amp, offset

def initParams(im,mode='basic',frac=.2,hdr=None):
    """Initial guess for beta, phi, nmax
    mode:
        basic: beta is determined based on the size of the image and the frac arguement, phi is 0
        fit: a 2D Gaussian is fit to the image, parameters derived from fit
            Theta_max = fit Gaussian width
            Theta_min = PSF FWHM
            beta ~ sqrt((Theta_max * Theta_min))
            phi ~ fit Gaussian rotation angle
            nmax ~ (Theta_max / Theta_min) + 1
    frac: fraction of image to use as the initial beta (basic)
    hdr: FITS header dictionary with PSF size (fit)
    fitting borrowed from: http://wiki.scipy.org/Cookbook/FittingData
    returns: beta, phi, nmax
    """
    if mode.startswith('basic'):
        return [frac*im.shape[0],frac*im.shape[1]], 0., int((frac/np.max(im.shape))-1)
    elif mode.startswith('fit'):
        #fit a 2D Gaussian to the image
        params=initGaussian(im)
        errorfunction = lambda p: np.ravel(ellipticalGaussian2D(*p)(*np.indices(im.shape)) - im)
        p, success = optimize.leastsq(errorfunction, params)
        Theta_max=2.3548*np.array([p[2],p[3]]) #FWHM

        #compute PSF size in pixels
        bpa=np.pi*hdr['bpa']/180.
        bmaj=np.pi*hdr['bmaj']/180.
        bmin=np.pi*hdr['bmin']/180.
        dra=np.pi*hdr['dra']/180.
        ddec=np.pi*hdr['ddec']/180.
        rotM=np.matrix([[np.cos(bpa),-1.*np.sin(bpa)],[np.sin(bpa),np.cos(bpa)]])
        rotDeltas=np.dot(rotM,np.array([dra,ddec])) #rotate delta RA and delta Dec
        psfPix=(np.array([bmaj,bmin])/rotDeltas)/2.3548 #go from FWHM to sigma, it is better to error on the side of higher order structure, then to miss it
        Theta_min=np.abs(np.array(psfPix).flatten())

        beta=np.sqrt(Theta_max*Theta_min)
        phi=p[4]
        print [(Theta_max[0]/Theta_min[0])-1,(Theta_max[1]/Theta_min[1])-1]
        nmax=[int((Theta_max[0]/Theta_min[0])+1),int((Theta_max[1]/Theta_min[1])+1)]
        return beta, phi, nmax

def initBetaPhi(im,mode='basic',frac=.2,circular=False):
    """Depreciated: use initParams
    Initial starting point for beta and phi
    mode:
        basic: beta is determined based on the size of the image and the frac arguement, phi is 0
        fit: a 2D Gaussian is fit to the image, parameters derived from fit
    frac: fraction of image to use as the initial beta (basic)
    circular: if True, return a circular fit parameter using the smallest width
    fitting borrowed from: http://wiki.scipy.org/Cookbook/FittingData
    """
    if mode.startswith('basic'):
        return [frac*im.shape[0],frac*im.shape[1]], 0.
    elif mode.startswith('fit'):
        params=initGaussian(im)
        errorfunction = lambda p: np.ravel(ellipticalGaussian2D(*p)(*np.indices(im.shape)) - im)
        p, success = optimize.leastsq(errorfunction, params)
        if circular:
            width=np.min([p[2],p[3]])
            return [width,width],p[4]
        else: return [p[2],p[3]],p[4]

def genPolarBasisFuncs(beta,nmax,phi,fourier=False):
    """Generate a list of basis functions
    beta: characteristic size of the shapelets (b_major, b_minor)
    nmax: maximum decomposition order
    phi: rotation angle
    fourier: return a FOurer transformed version of the basis functions
    """
    bfs=[]
    if type(nmax) is int: nmax=[nmax,nmax]
    for nn in range(nmax[0]):
        for mm in np.arange(-1*nn,nn+1):
            if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                bfs.append(shapelet.polarDimBasis(nn,mm,beta=beta,phi=phi,fourier=fourier))
    return bfs

def genPolarBasisMatrix(beta,nmax,phi,r,th,fourier=False):
    """Generate the n x k matrix of basis functions(k) for each pixel(n)
    beta: characteristic size of the shapelets (b_major, b_minor)
    nmax: maximum decomposition order
    phi: rotation angle
    r: radius matrix of n pixels
    th: theta matrix of n pixels
    fourier: return a FOurer transformed version of the basis functions
    """
    bvals=[]
    if type(nmax) is int: nmax=[nmax,nmax]
    for nn in range(nmax[0]):
        for mm in np.arange(-1*nn,nn+1):
            if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                bf=shapelet.polarDimBasis(nn,mm,beta=beta,phi=phi,fourier=fourier)
                bvals.append(shapelet.computeBasisPolar(bf,r,th).flatten())
    bm=np.array(bvals)
    return bm.transpose()

def genBasisFuncs(beta,nmax,phi,fourier=False):
    """Generate a list of basis functions
    beta: characteristic size of the shapelets (b_major, b_minor)
    nmax: maximum decomposition order
    phi: rotation angle
    fourier: return a FOurer transformed version of the basis functions
    """
    bfs=[]
    if type(nmax) is int: nmax=[nmax,nmax]
    for x in range(nmax[0]):
        for y in range(nmax[1]):
            bfs.append(shapelet.dimBasis2d(x,y,beta=beta,phi=phi,fourier=fourier))
    return bfs

def genBasisMatrix(beta,nmax,phi,rx,ry,fourier=False):
    """Generate the n x k matrix of basis functions(k) for each pixel(n)
    nmax: maximum decompisition order
    beta: characteristic size of the shapelet
    phi: rotation angle
    rx: range of x values to evaluate basis functions
    ry: range of y values to evaluate basis functions
    fourier: return a FOurer transformed version of the basis functions
    """
    bvals=[]
    if type(nmax) is int: nmax=[nmax,nmax]
    for x in range(nmax[0]):
        for y in range(nmax[1]):
            bf=shapelet.dimBasis2d(x,y,beta=beta,phi=phi,fourier=fourier)
            bvals.append(shapelet.computeBasis2d(bf,rx,ry).flatten())
    bm=np.array(bvals)
    return bm.transpose()

def solveCoeffs(m,im):
    """Solve for the basis function coefficients Theta^ using a Least Squares Esimation (Lupton pg. 84)
    theta^ = (m^T * m)^-1 * im
    n: number of pixels in the image
    k: number of basis functions (in all dimensions)
    m: n x k matrix of the basis functions for each pixel
    im: image of size n pixels
    returns theta_hat, a size k array of basis function coefficents"""
    #im_flat=np.reshape(im,im.shape[0]*im.shape[1],1)
    im_flat=im.flatten()
    mTm=np.dot(m.T,m)                    #matrix multiply m with it's transpose
    mTm_inv=np.linalg.inv(mTm)           #invert the k x k matrix
    mTm_inv_mT=np.dot(mTm_inv,m.T)       #matrix multiply the result with the transpose of m
    theta_hat=np.dot(mTm_inv_mT,im_flat) #compute the coefficents for the basis functions
    return theta_hat

def chi2PolarFunc(params,nmax,im,nm,order=['beta0','beta1','phi','xc','yc'],set_beta=[None,None],set_phi=None,set_xc=[None,None],r=None,th=None):
    """Function which is to be minimized in the chi^2 analysis for Polar shapelets
    params = [beta0, beta1, phi, xc, yc] or some subset
        beta0: characteristic size of shapelets, fit parameter
        beta1: characteristic size of shapelets, fit parameter
        phi: rotation angle of shapelets, fit parameter
        xc: x centroid of shapelets, fit parameter
        yc: y centroid of shapelets, fit parameter
    nmax: number of coefficents to use in the Laguerre polynomials
    im: observed image
    nm: noise map
    order: order of parameters
    fixed parameters: set_beta, set_phi, set_xc
    r: radius from centroid, array of im.shape, not required if xc and yc being fit
    th: angle from centroid, array of im.shape, not required if xc and yc being fit
    """
    #determine which parameters are being fit for, and which are not
    beta0=set_beta[0]
    beta1=set_beta[1]
    phi=set_phi
    xc=set_xc[0]
    yc=set_xc[1]
    fitParams={'beta':False,'phi':False,'xc':False}
    for pid,paramName in enumerate(order):
        if paramName=='beta0':
            beta0=params[pid]
            fitParams['beta']=True
        elif paramName=='beta1':
            beta1=params[pid]
            fitParams['beta']=True
        elif paramName=='phi':
            phi=params[pid]
            fitParams['phi']=True
        elif paramName=='xc':
            xc=params[pid]
            fitParams['xc']=True
        elif paramName=='yc':   
            yc=params[pid]
            fitParams['xc']=True

    if beta0<0.:
        print 'warning: beta going negative, setting to 0.0'
        beta0=0.
    if beta1<0.:
        print 'warning: beta going negative, setting to 0.0'
        beta1=0.
    print 'beta: (%f,%f)\t phi: %f\txc: (%f,%f)'%(beta0,beta1,phi,xc,yc)

    size=im.shape
    if fitParams['xc'] or r is None:
        r,th=shapelet.polarArray([xc,yc],size) #the redius,theta pairs need to updated if fitting for the xc centre or if not using the r,th inputs
    bvals=genPolarBasisMatrix([beta0,beta1],nmax,phi,r,th)
    coeffs=solveCoeffs(bvals,im)
    mdl=np.abs(img.constructModel(bvals,coeffs,size))
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2nmaxPolarFunc(params,im,nm,beta0,beta1,phi,xc):
    """
    params = [nmax]
        nmax: number of coefficents
    im: observed image
    nm: noise map
    beta0: fit beta value
    beta1: fit beta value
    phi: rotation angle
    xc: fit centroid position
    """
    #nmax=params[0]
    nmax=[params, params]
    size=im.shape
    r,th=shapelet.polarArray(xc,size)
    bvals=genPolarBasisMatrix([beta0,beta1],nmax,phi,r,th)
    coeffs=solveCoeffs(bvals,im)
    mdl=np.abs(img.constructModel(bvals,coeffs,size))
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2Func(params,nmax,im,nm,order=['beta0','beta1','phi','xc','yc'],set_beta=[None,None],set_phi=None,set_xc=[None,None],xx=None,yy=None):
    """Function which is to be minimized in the chi^2 analysis for Cartesian shapelets
    params = [beta0, beta1, phi, xc, yc] or some subset
        beta0: characteristic size of shapelets, fit parameter
        beta1: characteristic size of shapelets, fit parameter
        phi: rotation angle of shapelets, fit parameter
        xc: x centroid of shapelets, fit parameter
        yc: y centroid of shapelets, fit parameter
    nmax: number of coefficents to use in the Hermite polynomials
    im: observed image
    nm: noise map
    order: order of parameters
    fixed parameters: set_beta, set_phi, set_xc
    xx: X postiion grid, array of im.shape, not required if xc and yc being fit
    yy: Y position grid, array of im.shape, not required if xc and yc being fit
    """
    #determine which parameters are being fit for, and which are not
    betaX=set_beta[0]
    betaY=set_beta[1]
    phi=set_phi
    xc=set_xc[0]
    yc=set_xc[1]
    fitParams={'beta':False,'phi':False,'xc':False}
    for pid,paramName in enumerate(order):
        if paramName=='beta0':
            betaX=params[pid]
            fitParams['beta']=True
        elif paramName=='beta1':
            betaY=params[pid]
            fitParams['beta']=True
        elif paramName=='phi':
            phi=params[pid]
            fitParams['phi']=True
        elif paramName=='xc':
            xc=params[pid]
            fitParams['xc']=True
        elif paramName=='yc':   
            yc=params[pid]
            fitParams['xc']=True

    if betaX<0.:
        print 'warning: beta going negative, setting to 0.0'
        betaX=0.
    if betaY<0.:
        print 'warning: beta going negative, setting to 0.0'
        betaY=0.
    print 'beta: (%f,%f)\t phi: %f\txc: (%f,%f)'%(betaX,betaY,phi,xc,yc)

    size=im.shape
    if fitParams['xc'] or xx is None:
        #shift the (0,0) point to the centroid
        rx=np.array(range(0,size[0]),dtype=float)-xc
        ry=np.array(range(0,size[1]),dtype=float)-yc
        xx,yy=shapelet.xy2Grid(rx,ry)
    bvals=genBasisMatrix([betaX,betaY],nmax,phi,xx,yy)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

def chi2nmaxFunc(params,im,nm,beta0,beta1,phi,xc):
    """
    params = [nmaxX,nmaxY]
        nmax: number of coefficents to use in x,y
    im: observed image
    nm: noise map
    beta0: characteristic size of shapelet
    beta1: characteristic size of shapelet
    phi: rotation angle of shapelets
    xc: fit centroid position
    """
    nmax=params
    size=im.shape
    #shift the (0,0) point to the centroid
    rx=np.array(range(0,size[0]),dtype=float)-xc[0]
    ry=np.array(range(0,size[1]),dtype=float)-xc[1]
    xx,yy=shapelet.xy2Grid(rx,ry)

    bvals=genBasisMatrix([beta0,beta1],[nmax,nmax],phi,xx,yy)
    coeffs=solveCoeffs(bvals,im)
    mdl=img.constructModel(bvals,coeffs,size)
    return np.sum((im-mdl)**2 / nm**2)/(size[0]*size[1])

if __name__ == "__main__":

    print '============================================'
    print 'Testing decomp module:'
    print '============================================'
    import fileio
    import sys
    write_files = False #Flag to write coeff files
    tc=0
    te=0
    
    im,hdr=fileio.readFITS('../data/N6251_test.fits',hdr=True)
    subim=img.selPxRange(im,[1028,1097,1025,1074])
    
    #guess beta
    #initBeta(im):
    tc+=1
    try:
        beta0,phi0=initBetaPhi(subim,mode='basic')
        print beta0, phi0
        beta0,phi0=initBetaPhi(subim,mode='fit')
        print beta0, phi0
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #genPolarBasisMatrix(beta,nmax,r,th):
    #solveCoeffs(m,im):
    tc+=1
    try:
        beta0,phi0=initBetaPhi(subim,mode='fit')
        xc=img.maxPos(subim)
        r0,th0=shapelet.polarArray(xc,subim.shape)
        mPolar=genPolarBasisMatrix(beta0,5,0.,r0,th0)
        coeffs=solveCoeffs(mPolar,subim)
        print coeffs
        if write_files: fileio.writeLageurreCoeffs('testLageurre.pkl',coeffs,xc,subim.shape,beta0,0.,[5,5],pos=[hdr['ra'],hdr['dec'],hdr['dra'],hdr['ddec']],info='Test Lageurre coeff file')
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #genBasisMatrix(beta,nmax,rx,ry):
    #solveCoeffs(m,im):
    tc+=1
    try:
        beta0,phi0=initBetaPhi(subim,mode='fit')
        xc=img.maxPos(subim)
        rx=np.array(range(0,subim.shape[0]),dtype=float)-xc[0]
        ry=np.array(range(0,subim.shape[1]),dtype=float)-xc[1]
        xx,yy=shapelet.xy2Grid(rx,ry)
        mCart=genBasisMatrix(beta0,[5,5],phi0,xx,yy)
        coeffs=solveCoeffs(mCart,subim)
        print coeffs
        if write_files: fileio.writeHermiteCoeffs('testHermite.pkl',coeffs,xc,subim.shape,beta0,0.,5,pos=[hdr['ra'],hdr['dec'],hdr['dra'],hdr['ddec']],info='Test Hermite coeff file')
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    beta0,phi0=initBetaPhi(subim,mode='fit')
    xc=img.maxPos(subim)
    mean,std=img.estimateNoise(subim,mode='basic')
    nm=img.makeNoiseMap(subim.shape,mean,std)
    nmax=[10,10]

    #chi2PolarFunc(params,nmax,im,nm):
    tc+=1
    try:
        print chi2PolarFunc([beta0[0],beta0[1],phi0,xc[0],xc[1]],nmax,subim,nm,order=['beta0','beta1','phi','xc','yc']) #fit: all
        print chi2PolarFunc([beta0[0],beta0[1]],nmax,subim,nm,order=['beta0','beta1'],set_phi=phi0,set_xc=xc) #fit: beta
        print chi2PolarFunc([phi0],nmax,subim,nm,order=['phi'],set_beta=beta0,set_xc=xc) #fit: phi
        print chi2PolarFunc([xc[0],xc[1]],nmax,subim,nm,order=['xc','yc'],set_beta=beta0,set_phi=phi0) #fit: xc
        print chi2PolarFunc([beta0[0],beta0[1],phi0],nmax,subim,nm,order=['beta0','beta1','phi'],set_xc=xc) #fit: beta, phi
        print chi2PolarFunc([beta0[0],beta0[1],xc[0],xc[1]],nmax,subim,nm,order=['beta0','beta1','xc','yc'],set_phi=phi0) #fit: beta, xc
        print chi2PolarFunc([phi0,xc[0],xc[1]],nmax,subim,nm,order=['phi','xc','yc'],set_beta=beta0) #fit: phi, xc
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #chi2nmaxPolarFunc(params,im,nm,beta,xc):
    tc+=1
    try:
        print chi2nmaxPolarFunc(5,subim,nm,beta0[0],beta0[1],0.,xc)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #chi2Func(params,nmax,im,nm):
    tc+=1
    try:
        print chi2Func([beta0[0],beta0[1],phi0,xc[0],xc[1]],nmax,subim,nm,order=['beta0','beta1','phi','xc','yc']) #fit: all
        print chi2Func([beta0[0],beta0[1]],nmax,subim,nm,order=['beta0','beta1'],set_phi=phi0,set_xc=xc) #fit: beta
        print chi2Func([phi0],nmax,subim,nm,order=['phi'],set_beta=beta0,set_xc=xc) #fit: phi
        print chi2Func([xc[0],xc[1]],nmax,subim,nm,order=['xc','yc'],set_beta=beta0,set_phi=phi0) #fit: xc
        print chi2Func([beta0[0],beta0[1],phi0],nmax,subim,nm,order=['beta0','beta1','phi'],set_xc=xc) #fit: beta, phi
        print chi2Func([beta0[0],beta0[1],xc[0],xc[1]],nmax,subim,nm,order=['beta0','beta1','xc','yc'],set_phi=phi0) #fit: beta, xc
        print chi2Func([phi0,xc[0],xc[1]],nmax,subim,nm,order=['phi','xc','yc'],set_beta=beta0) #fit: phi, xc
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #chi2nmaxFunc(params,im,nm,beta,xc):
    tc+=1
    try:
        print chi2nmaxFunc(5,subim,nm,beta0[0],beta0[1],0.,xc)
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    print '============================================'
    print '%i of %i tests succeeded'%(tc-te,tc)
    print '============================================'

