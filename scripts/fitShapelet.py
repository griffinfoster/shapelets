#!/usr/bin/env python
"""
Fit for shapelet coefficients across beta, xc, phi, and n_max
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches
from scipy import optimize
import shapelets

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_IMAGE')
    o.set_description(__doc__)
    o.add_option('-r', '--region', dest='region', default=None,
        help='Region of image to decompose into shapelets, (xmin,xmax,ymin,ymax), default: None')
    o.add_option('-N', '--noise_region', dest='nregion', default=None,
        help='Region of image to use to create a noise map, if set to None the entire image is used, this is not used in the script, (xmin,xmax,ymin,ymax), default: None')
    o.add_option('-m', '--mode', dest='mode', default='cart',
        help='Set the shapelet mode, cartesian or polar, default: cartesian')
    o.add_option('-o', '--outfile', dest='ofn', default='shapeletCoeffs.pkl',
        help='Coefficients output filename, default: shapeletCoeffs.pkl')
    o.add_option('--max', dest='max_pos', action="store_true", default=False,
        help='Override centroid position to be the position of max intensity')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save the figure, requires filename')

    o.add_option('-x','--init_xc', dest='init_xc', default=None,
        help='Initial parameter: set a x,y pixel position for initial center, if using a region it is based on the relative position, default: centroid of image/region')
    o.add_option('--set_xc', dest='set_xc', action='store_true',
        help='Set parameter: set init_xc x,y pixel position for center, these parameters will not be fit for if set')
    o.add_option('-b', '--init_beta', dest='init_beta', default=None,
        help='Initial parameter: initial beta value, can be two values i.e. \'25.0,30.5\', default: None, guess is made based on Gaussian fit')
    o.add_option('--set_beta', dest='set_beta', action='store_true',
        help='Set parameter: set init_beta beta value, these parameters will not be fit for if set')
    o.add_option('-p','--init_phi', dest='init_phi', default=None,
        help='Initial parameter: inital rotation angle (radians), only used when beta is manually input, default: 0')
    o.add_option('--set_phi', dest='set_phi', action='store_true',
        help='Set parameter: set init_phi rotation angle, this parameter will not be fit for if set')
    o.add_option('-n', '--nmax', dest='nmax', default='5',
        help='Size of coefficient dimensions for minimization fit, can be two values i.e. \'4,5\', default: 5')

    o.add_option('--fitter',dest='fitterMethod', default='Nelder-Mead',
        help='Fitting method: Nelder-Mead, Powell, CG, BFGS, see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.minimize.html, default: Nelder-Mead')
    o.add_option('--xtol', dest='xtol', default=0.001, type='float',
        help='Relative error in parameters acceptable for convergence, default: 0.001')
    o.add_option('--ftol', dest='ftol', default=0.001, type='float',
        help='Relative error in chi^2 function acceptable for convergence, default: 0.001')
    o.add_option('--maxiter', dest='maxiter', default=250, type='int',
        help='Maximum number of iterations to perform, default: 250')
    o.add_option('-B', '--brute', dest='brute', default=15, type='int',
        help='Maximum basis function order to use when running brute force method, default: 15')
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    im0=shapelets.fileio.readFITS(ifn)
    extent=[0,im0.shape[0],0,im0.shape[1]]
    if not (opts.region is None):
        extent=map(int, opts.region.split(','))
        im=shapelets.img.selPxRange(im0,extent)
    else:
        im=im

    #noise map
    if opts.nregion is None:
        #use the image region for noise estimation
        mean,std=shapelets.img.estimateNoise(im,mode='basic')
        nm=shapelets.img.makeNoiseMap(im.shape,mean,std)
    else:
        #use a specific region for noise estimation
        nextent=map(int, opts.nregion.split(','))
        mean,std=shapelets.img.estimateNoise(shapelets.img.selPxRange(im0,nextent),mode='basic')
        nm=shapelets.img.makeNoiseMap(im.shape,mean,std)

    #determine set parameters
    set_xc=opts.set_xc
    set_beta=opts.set_beta
    set_phi=opts.set_phi

    #select initial beta, phi, and xc
    if opts.init_beta==None:
        beta0,phiTemp=shapelets.decomp.initBetaPhi(im,mode='fit')
    else:
        beta0=map(float,opts.init_beta.split(','))
        if len(beta0)==1:
            beta0=[beta0[0],beta0[0]]
        else:
            beta0=[beta0[0],beta0[1]]

    if opts.init_phi==None:
        betaTemp,phi0=shapelets.decomp.initBetaPhi(im,mode='fit')
    else:
        phi0=float(opts.init_phi)
    
    if opts.init_xc==None:
        xc=shapelets.img.centroid(im)
    else:
        xc=map(float,opts.init_xc.split(','))
        #correct for position if using only a region of the image
        #xc[0]-=extent[0]
        #xc[1]-=extent[2]

    nmax=opts.nmax.split(',')
    if len(nmax)==1:
        nmax=[int(nmax[0])+1,int(nmax[0])+1]
    else:
        nmax=[int(nmax[0])+1,int(nmax[1])+1]

    print 'Using beta: (%f,%f) :: \tphi: %f radians :: \tcentre: x,y=(%f,%f) :: \tnmax: (%i,%i)'%(beta0[0],beta0[1],phi0,xc[0],xc[1],nmax[0]-1,nmax[1]-1)
    print 'Fitting xc : %r\nFitting beta : %r\nFitting phi : %r'%(not(set_xc),not(set_beta),not(set_phi))

    if opts.mode.startswith('pol'):
        r0,th0=shapelets.shapelet.polarArray(xc,im.shape)

        #scipy-based minimizer
        if set_xc:
            if set_beta:
                if set_phi:
                    #same as solveShapelets, no minimization
                    print 'No parameters to minimize, solving for coefficients with input values'
                    beta1=beta0
                    phi1=phi0
                    xc1=xc
                else:
                    print 'Running minimization for phi only...'
                    res=optimize.minimize(shapelets.decomp.chi2PolarFunc,[phi0],args=(nmax,im,nm,['phi'],beta0,None,xc,r0,th0),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=beta0
                    phi1=res['x'][0]
                    xc1=xc
            else:
                if set_phi:
                    print 'Running minimization for beta only...'
                    res=optimize.minimize(shapelets.decomp.chi2PolarFunc,[beta0[0], beta0[1]],args=(nmax,im,nm,['beta0','beta1'],[None,None],phi0,xc,r0,th0),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=[res['x'][0],res['x'][1]]
                    phi1=phi0
                    xc1=xc
                else:
                    print 'Running minimization for beta and phi...'
                    res=optimize.minimize(shapelets.decomp.chi2PolarFunc,[beta0[0], beta0[1], phi0],args=(nmax,im,nm,['beta0','beta1','phi'],[None,None],None,xc,r0,th0),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=[res['x'][0],res['x'][1]]
                    phi1=res['x'][2]
                    xc1=xc
        else:
            if set_beta:
                if set_phi:
                    print 'Running minimization for centroid only...'
                    res=optimize.minimize(shapelets.decomp.chi2PolarFunc,[xc[0],xc[1]],args=(nmax,im,nm,['xc','yc'],beta0,phi0,[None,None],None,None),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=beta0
                    phi1=phi0
                    xc1=[res['x'][0],res['x'][1]]
                else:
                    print 'Running minimization for phi and centroid...'
                    res=optimize.minimize(shapelets.decomp.chi2PolarFunc,[phi0,xc[0],xc[1]],args=(nmax,im,nm,['phi','xc','yc'],beta0,None,[None,None],None,None),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=beta0
                    phi1=res['x'][0]
                    xc1=[res['x'][1],res['x'][2]]
            else:
                if set_phi:
                    print 'Running minimization for beta and centroid...'
                    res=optimize.minimize(shapelets.decomp.chi2PolarFunc,[beta0[0],beta0[1],xc[0],xc[1]],args=(nmax,im,nm,['beta0','beta1','xc','yc'],[None,None],phi0,[None,None],None,None),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=[res['x'][0],res['x'][1]]
                    phi1=phi0
                    xc1=[res['x'][2],res['x'][3]]
                else:
                    print 'Running minimization for beta, phi and centroid...'
                    res=optimize.minimize(shapelets.decomp.chi2PolarFunc,[beta0[0], beta0[1], phi0, xc[0], xc[1]],args=(nmax,im,nm),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=[res['x'][0],res['x'][1]]
                    phi1=res['x'][2]
                    xc1=[res['x'][3],res['x'][4]]
        print '\tDone'

        #scipy optimize brute force over a range of N values
        n0=1
        n1=opts.brute+1
        print 'Running brute force for size of N on range [%i:%i]...'%(n0,n1-1)
        nmax1=optimize.brute(shapelets.decomp.chi2nmaxPolarFunc,[np.s_[n0:n1:1]],args=(im,nm,beta1[0],beta1[1],phi1,xc1),finish=None)
        nmax1=[int(nmax1),int(nmax1)]
        print 'Using %i x %i coefficients'%(nmax1[0],nmax1[1])
        print '\tDone'

        print 'Solution:'
        print '\tbeta: (%f,%f) \tphi: %f rad \tcentroid: (%f,%f) pixels \t ncoeffs: %i x %i'%(beta1[0], beta1[1], phi1, xc1[0], xc1[1], nmax1[0],nmax1[1])

        #plot: data, model, residual: model-data, coeffs
        fig = plt.figure()
        ax = fig.add_subplot(221)
        plt.title('Image')
        plt.imshow(im)
        e=matplotlib.patches.Ellipse(xy=xc,width=2.*np.max(beta0),height=2.*np.min(beta0),angle=(180.*phi0/np.pi))
        e.set_clip_box(ax.bbox)
        e.set_alpha(0.3)
        e.set_facecolor('black')
        ax.add_artist(e)
        plt.text(xc[0],xc[1],'+',horizontalalignment='center',verticalalignment='center')
        plt.colorbar()
        
        plt.subplot(222)
        plt.title('Model')
        r1,th1=shapelets.shapelet.polarArray(xc1,im.shape)
        bvals=shapelets.decomp.genPolarBasisMatrix(beta1,nmax1,phi1,r1,th1)
        coeffs=shapelets.decomp.solveCoeffs(bvals,im)
        mdl=np.abs(shapelets.img.constructModel(bvals,coeffs,im.shape))
        plt.imshow(mdl)
        plt.colorbar()
        
        plt.subplot(223)
        plt.title('Residual')
        res=im-mdl
        plt.imshow(res)
        plt.colorbar()
        
        plt.subplot(224)
        plt.title('Coefficients')
        cimR=shapelets.img.polarCoeffImg(coeffs.real,nmax1)
        cimI=shapelets.img.polarCoeffImg(coeffs.imag,nmax1)
        cimI=np.fliplr(cimI)
        cim=np.concatenate((cimR,cimI),axis=1)
        plt.pcolor(cim)
        plt.colorbar()

        ofn=opts.ofn
        print 'Writing to file:',ofn
        shapelets.fileio.writeLageurreCoeffs(ofn,coeffs,xc,im.shape,beta0,phi0,nmax1,info=ifn)
        
    else:
        rx=np.array(range(0,im.shape[0]),dtype=float)-xc[0]
        ry=np.array(range(0,im.shape[1]),dtype=float)-xc[1]
        xx0,yy0=shapelets.shapelet.xy2Grid(rx,ry)

        #scipy-based minimizer
        if set_xc:
            if set_beta:
                if set_phi:
                    #same as solveShapelets, no minimization
                    print 'No parameters to minimize, solving for coefficients with input values'
                    beta1=beta0
                    phi1=phi0
                    xc1=xc
                else:
                    print 'Running minimization for phi only...'
                    res=optimize.minimize(shapelets.decomp.chi2Func,[phi0],args=(nmax,im,nm,['phi'],beta0,None,xc,xx0,yy0),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=beta0
                    phi1=res['x'][0]
                    xc1=xc
            else:
                if set_phi:
                    print 'Running minimization for beta only...'
                    res=optimize.minimize(shapelets.decomp.chi2Func,[beta0[0], beta0[1]],args=(nmax,im,nm,['beta0','beta1'],[None,None],phi0,xc,xx0,yy0),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=[res['x'][0],res['x'][1]]
                    phi1=phi0
                    xc1=xc
                else:
                    print 'Running minimization for beta and phi...'
                    res=optimize.minimize(shapelets.decomp.chi2Func,[beta0[0], beta0[1], phi0],args=(nmax,im,nm,['beta0','beta1','phi'],[None,None],None,xc,xx0,yy0),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=[res['x'][0],res['x'][1]]
                    phi1=res['x'][2]
                    xc1=xc
        else:
            if set_beta:
                if set_phi:
                    print 'Running minimization for centroid only...'
                    res=optimize.minimize(shapelets.decomp.chi2Func,[xc[0],xc[1]],args=(nmax,im,nm,['xc','yc'],beta0,phi0,[None,None],None,None),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=beta0
                    phi1=phi0
                    xc1=[res['x'][0],res['x'][1]]
                else:
                    print 'Running minimization for phi and centroid...'
                    res=optimize.minimize(shapelets.decomp.chi2Func,[phi0,xc[0],xc[1]],args=(nmax,im,nm,['phi','xc','yc'],beta0,None,[None,None],None,None),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=beta0
                    phi1=res['x'][0]
                    xc1=[res['x'][1],res['x'][2]]
            else:
                if set_phi:
                    print 'Running minimization for beta and centroid...'
                    res=optimize.minimize(shapelets.decomp.chi2Func,[beta0[0],beta0[1],xc[0],xc[1]],args=(nmax,im,nm,['beta0','beta1','xc','yc'],[None,None],phi0,[None,None],None,None),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=[res['x'][0],res['x'][1]]
                    phi1=phi0
                    xc1=[res['x'][2],res['x'][3]]
                else:
                    print 'Running minimization for beta, phi and centroid...'
                    res=optimize.minimize(shapelets.decomp.chi2Func,[beta0[0], beta0[1], phi0, xc[0], xc[1]],args=(nmax,im,nm),method=opts.fitterMethod,options={'xtol':opts.xtol,'ftol':opts.ftol,'maxiter':opts.maxiter})
                    print res
                    beta1=[res['x'][0],res['x'][1]]
                    phi1=res['x'][2]
                    xc1=[res['x'][3],res['x'][4]]
        print '\tDone'

        #scipy optimize brute force over a range of N values
        n0=1
        n1=opts.brute+1
        print 'Running brute force for size of N on range [%i:%i]...'%(n0,n1-1)
        nmax1=optimize.brute(shapelets.decomp.chi2nmaxFunc,[np.s_[n0:n1:1]],args=(im,nm,beta1[0],beta1[1],phi1,xc1),finish=None)
        nmax1=[int(nmax1),int(nmax1)]
        print 'Using %i x %i coefficients'%(nmax1[0],nmax1[1])
        print '\tDone'

        print 'Solution:'
        print '\tbeta: (%f,%f) \tphi: %f rad \tcentroid: (%f,%f) pixels \t ncoeffs: %i x %i'%(beta1[0], beta1[1], phi1, xc1[0], xc1[1], nmax1[0],nmax1[1])

        #plot: data, model, residual: model-data, coeffs
        fig = plt.figure()
        ax = fig.add_subplot(221)
        plt.title('Image')
        plt.imshow(im)
        e=matplotlib.patches.Ellipse(xy=xc,width=2.*np.max(beta0),height=2.*np.min(beta0),angle=(180.*phi0/np.pi))
        e.set_clip_box(ax.bbox)
        e.set_alpha(0.3)
        e.set_facecolor('black')
        ax.add_artist(e)
        plt.text(xc[0],xc[1],'+',horizontalalignment='center',verticalalignment='center')
        plt.colorbar()
        
        plt.subplot(222)
        plt.title('Model')
        rx=np.array(range(0,im.shape[0]),dtype=float)-xc1[0]
        ry=np.array(range(0,im.shape[1]),dtype=float)-xc1[1]
        xx,yy=shapelets.shapelet.xy2Grid(rx,ry)
        bvals=shapelets.decomp.genBasisMatrix(beta1,nmax1,phi1,xx,yy)
        coeffs=shapelets.decomp.solveCoeffs(bvals,im)
        mdl=shapelets.img.constructModel(bvals,coeffs,im.shape)
        plt.imshow(mdl)
        plt.text(xc[0],xc[1],'+',horizontalalignment='center',verticalalignment='center')
        plt.colorbar()
        
        plt.subplot(223)
        plt.title('Residual')
        res=im-mdl
        plt.imshow(res)
        plt.colorbar()

        plt.subplot(224)
        plt.title('Coefficients')
        sqCoeffs=np.reshape(coeffs,nmax1)
        plt.pcolor(sqCoeffs)
        plt.colorbar()
        
        ofn=opts.ofn
        print 'Writing to file:',ofn
        shapelets.fileio.writeHermiteCoeffs(ofn,coeffs,xc1,im.shape,beta1,phi1,nmax1,info=ifn)
        
    if not (opts.savefig is None):
        plt.savefig(opts.savefig)
    else: plt.show()
