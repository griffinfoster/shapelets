#!/usr/bin/env python
"""
Plot a Hermite shapelet coefficient file

ToDo:
    Fourier Transform for Polar Shapelets
    Understand size,res rescaling in the Fourier Transform
"""

import sys,os
import numpy as n
import pylab as p
import shapelets

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] COEFF_FILE')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    ifn=args[0]
    d=shapelets.fileio.readCoeffs(ifn)

    if d['mode'].startswith('herm'):
        #model
        rx=n.array(range(0,d['size'][0]),dtype=float)-d['xc'][0]
        ry=n.array(range(0,d['size'][1]),dtype=float)-d['xc'][1]
        #print rx
        #print ry
        #print len(rx),len(ry)
        bvals=shapelets.decomp.genBasisMatrix(d['beta'],d['norder'],rx,ry)
        mdl=shapelets.img.constructModel(bvals,d['coeffs'],d['xc'],d['size'])
        
        #ft
        usize=(float(len(ry))/float(len(rx)))*1./n.abs(rx[0]-rx[1]) #!!!
        vsize=(float(len(rx))/float(len(ry)))*1./n.abs(ry[0]-ry[1]) #!!!
        #print usize,vsize
        ures=(float(len(ry))/float(len(rx)))*1./n.abs(rx[0]-rx[-1]) #!!!
        vres=(float(len(rx))/float(len(ry)))*1./n.abs(ry[0]-ry[-1]) #!!!
        #print ures,vres
        umin=ures*(rx[0]+.5)
        vmin=vres*(ry[0]+.5)
        uu=n.arange(0,usize+ures,ures)+umin
        vv=n.arange(0,vsize+vres,vres)+vmin

        usize=5.*1./n.abs(rx[0]-rx[1]) #!!!
        vsize=5.*1./n.abs(ry[0]-ry[1]) #!!!
        ures=5.*1./n.abs(rx[0]-rx[-1]) #!!!
        vres=5.*1./n.abs(ry[0]-ry[-1]) #!!!
        uu=n.arange(-.5*usize,.5*usize,ures)
        vv=n.arange(-.5*vsize,.5*vsize,vres)
        #print len(uu),len(vv)
        #print uu,vv
        bvals=[]
        basis=shapelets.shapelet.ftHermiteBasis(d['beta'],d['norder'])
        for bf in basis:
            bvals.append(shapelets.shapelet.computeBasis2d(bf,uu,vv).flatten())
        bm=n.array(bvals)
        mdl_ft=shapelets.img.constructModel(bm.transpose(),d['coeffs'],d['xc'],[len(uu),len(vv)])
        #mdl_ft=mdl_ft.real

        #uu=n.array(range(0,d['size'][0]),dtype=float)-d['xc'][0]
        #vv=n.array(range(0,d['size'][1]),dtype=float)-d['xc'][1]
        ##print uu
        ##print vv
        #bvals=[]
        #basis=shapelets.shapelet.ftHermiteBasis([1./d['beta'][0],1./d['beta'][1]],d['norder'])
        #for bf in basis:
        #    bvals.append(shapelets.shapelet.computeBasis2d(bf,uu,vv).flatten())
        #bm=n.array(bvals)
        #mdl_ft=shapelets.img.constructModel(bm.transpose(),d['coeffs'],[len(uu),len(vv)])
        #mdl_ft=mdl_ft.real

        p.suptitle('Hermite')
    elif d['mode'].startswith('lag'):
        r0,th0=shapelets.shapelet.polarArray(d['xc'],d['size'])
        bvals=shapelets.decomp.genPolarBasisMatrix(d['beta'],d['norder'],r0,th0)
        mdl=n.abs(shapelets.img.constructModel(bvals,d['coeffs'],d['size']))
        coeffs=d['coeffs']
        cimR=shapelets.img.polarCoeffImg(coeffs.real,d['norder'])
        cimI=shapelets.img.polarCoeffImg(coeffs.imag,d['norder'])
        cimI=n.fliplr(cimI)
        cim=n.concatenate((cimR,cimI),axis=1)
        p.suptitle('Laguerre')
    
    p.subplot(231)
    p.title('Model')
    p.imshow(mdl)
    p.text(d['xc'][1],d['xc'][0],'+')
    #p.colorbar()

    p.subplot(234)
    mdl_fft2=n.fft.fft2(mdl)
    axis0_shift=int(mdl_fft2.shape[0]/2)
    axis1_shift=int(mdl_fft2.shape[1]/2)
    mdl_fft2=n.concatenate((mdl_fft2[axis0_shift:,:],mdl_fft2[:axis0_shift,:]),axis=0)
    mdl_fft2=n.concatenate((mdl_fft2[:,axis1_shift:],mdl_fft2[:,:axis1_shift]),axis=1)
    p.imshow(n.abs(mdl_fft2))

    p.subplot(232)
    p.title('Fourier Transform')
    #p.imshow(mdl_ft.real)
    p.imshow(n.abs(mdl_ft))
    #p.colorbar()

    p.subplot(235)
    mdl_ft_fft2=n.fft.fft2(mdl_ft)
    axis0_shift=int(mdl_ft_fft2.shape[0]/2)
    axis1_shift=int(mdl_ft_fft2.shape[1]/2)
    mdl_ft_fft2=n.concatenate((mdl_ft_fft2[axis0_shift:,:],mdl_ft_fft2[:axis0_shift,:]),axis=0)
    mdl_ft_fft2=n.concatenate((mdl_ft_fft2[:,axis1_shift:],mdl_ft_fft2[:,:axis1_shift]),axis=1)
    p.imshow(n.abs(mdl_ft_fft2))

    p.subplot(233)
    p.title('Coefficents')
    if d['mode'].startswith('herm'):
        coeffs=n.reshape(d['coeffs'],d['norder'])
        p.pcolor(coeffs)
    elif d['mode'].startswith('lag'):
        p.pcolor(cim)
    p.colorbar()
    
    p.show()
