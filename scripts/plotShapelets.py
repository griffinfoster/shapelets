#!/usr/bin/env python
"""
Plot a set of shapelet basis functions
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
import shapelets

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options]')
    o.set_description(__doc__)
    o.add_option('-n', '--nmax', dest='nmax', default='5',
        help='Maximum order of shapelet decomposition basis functions, if using Hermite can use 2 values i.e. \'4,5\', default: 5')
    o.add_option('-p', '--polar', dest='polar', action='store_true',
        help='Put in polar shapelet mode, default: False (cartesian)')
    o.add_option('-b', '--beta', dest='beta', default='1.0',
        help='Characteristic shapelet size, can use 1(symetric beta), 2(elliptical beta), 3(elliptical beta with rotation) values i.e. \'2.1,3.5\' default: 1.0')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save the figure, requires filename')
    o.add_option('-f', '--ft', dest='fourier', default=None,
        help='Perform a Fourier transform on the basis functions, mode: fft (compute 2D Fourer transform). scale (compute the Fourier transform as a rescaling)')
    opts, args = o.parse_args(sys.argv[1:])

    nmax=opts.nmax.split(',')
    if len(nmax)==1: nmax=[int(nmax[0])+1,int(nmax[0])+1]
    else: nmax=[int(nmax[0])+1,int(nmax[1])+1]

    betaphi=opts.beta.split(',')
    if len(betaphi)==1:
        beta=[float(betaphi[0]),float(betaphi[0])]
        phi=0.
    elif len(betaphi)==2:
        beta=[float(betaphi[0]),float(betaphi[1])]
        phi=0.
    else:
        beta=[float(betaphi[0]),float(betaphi[1])]
        phi=float(betaphi[2])

    xlim=[-5,5]
    ylim=[-5,5]
    rx=np.linspace(xlim[0],xlim[1],num=128)
    ry=np.linspace(ylim[0],ylim[1],num=128)

    if opts.polar:
        nmax=nmax[0]
        print 'polar shapelets'
        fullImgReal=np.zeros((len(ry)*nmax*2,len(rx)*nmax))
        fullImgImag=np.zeros((len(ry)*nmax*2,len(rx)*nmax))
        yOffset=len(ry)*nmax
        r,th=shapelets.shapelet.xy2rthGrid(rx,ry)
        for nn in range(nmax):
            for mm in np.arange(-1*nn,nn+1):
                if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                    if opts.fourier is None:
                        bf=shapelets.shapelet.polarDimBasis(nn,mm,beta=beta,phi=phi)
                        bval=shapelets.shapelet.computeBasisPolar(bf,r,th)
                    elif opts.fourier.startswith('fft'):
                        bf=shapelets.shapelet.polarDimBasis(nn,mm,beta=beta,phi=phi)
                        bval=shapelets.shapelet.computeBasisPolar(bf,r,th)
                        bval=np.fft.fftshift(np.fft.fft2(np.fft.fftshift(bval)))
                    elif opts.fourier.startswith('scale'):
                        bf=shapelets.shapelet.polarDimBasis(nn,mm,beta=beta,phi=phi,fourier=True)
                        bval=shapelets.shapelet.computeBasisPolar(bf,r,th)
                    fullImgReal[mm*len(ry)+yOffset:(mm+1)*len(ry)+yOffset,nn*len(rx):(nn+1)*len(rx)]=bval.real
                    fullImgImag[mm*len(ry)+yOffset:(mm+1)*len(ry)+yOffset,nn*len(rx):(nn+1)*len(rx)]=bval.imag
        fig=plt.figure()
        fig.subplots_adjust(wspace=0.3)
        
        plt.subplot(121)
        plt.imshow(fullImgReal)
        for nn in range(nmax):
            for mm in np.arange(-1*nn,nn+1):
                if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                    plt.fill([nn*len(rx),(nn+1)*len(rx),(nn+1)*len(rx),nn*len(rx)],[mm*len(ry)+yOffset,mm*len(ry)+yOffset,(mm+1)*len(ry)+yOffset,(mm+1)*len(ry)+yOffset],fill=False)
        plt.xlim(xmin=0,xmax=len(rx)*nmax)
        plt.ylim(ymin=len(ry)*nmax*2,ymax=len(ry))
        plt.xticks(np.arange(0,len(rx)*nmax,len(rx))+(len(rx)/2),range(nmax))
        plt.yticks(np.arange(0,len(ry)*(nmax*2+1),len(ry)+1)+(len(ry)/2),nmax-np.arange(nmax*2+1))
        plt.title('Real')
        plt.colorbar()
        
        plt.subplot(122)
        plt.imshow(fullImgImag)
        for nn in range(nmax):
            for mm in np.arange(-1*nn,nn+1):
                if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                    plt.fill([nn*len(rx),(nn+1)*len(rx),(nn+1)*len(rx),nn*len(rx)],[mm*len(ry)+yOffset,mm*len(ry)+yOffset,(mm+1)*len(ry)+yOffset,(mm+1)*len(ry)+yOffset],fill=False)
        plt.xlim(xmin=0,xmax=len(rx)*nmax)
        plt.ylim(ymin=len(ry)*nmax*2,ymax=len(ry))
        plt.xticks(np.arange(0,len(rx)*nmax,len(rx))+(len(rx)/2),range(nmax))
        plt.yticks(np.arange(0,len(ry)*(nmax*2+1),len(ry)+1)+(len(ry)/2),nmax-np.arange(nmax*2+1))
        plt.title('Imaginary')
        plt.colorbar()
       
        plt.suptitle('Lageurre Basis Functions (Polar)')
    
    else:
        print 'cartesian shapelets'
        fig=plt.figure()
        fullImg=np.zeros((len(rx)*nmax[1],len(ry)*nmax[0]))
        xx,yy=shapelets.shapelet.xy2Grid(rx,ry)
        for n0 in range(nmax[1]):
            for n1 in range(nmax[0]):
                if opts.fourier is None:
                    bf=shapelets.shapelet.dimBasis2d(n0,n1,beta=beta,phi=phi)
                    bval=shapelets.shapelet.computeBasis2d(bf,xx,yy)
                elif opts.fourier.startswith('fft'):
                    bf=shapelets.shapelet.dimBasis2d(n0,n1,beta=beta,phi=phi)
                    bval=shapelets.shapelet.computeBasis2d(bf,xx,yy)
                    bval=np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(bval))))
                elif opts.fourier.startswith('scale'):
                    bf=shapelets.shapelet.dimBasis2d(n0,n1,beta=beta,phi=phi,fourier=True)
                    bval=shapelets.shapelet.computeBasis2d(bf,xx,yy)
                fullImg[n0*len(rx):(n0+1)*len(rx),n1*len(ry):(n1+1)*len(ry)]=bval
        
        plt.imshow(fullImg)
        for n0 in range(nmax[1]):
            for n1 in range(nmax[0]):
                plt.fill([n1*len(rx),(n1+1)*len(rx),(n1+1)*len(rx),n1*len(rx)],[n0*len(ry),n0*len(ry),(n0+1)*len(ry),(n0+1)*len(ry)],fill=False)
        plt.xlim(xmin=0,xmax=len(rx)*nmax[0])
        plt.xticks(np.arange(0,len(rx)*nmax[0],len(rx))+(len(rx)/2),range(nmax[0]))
        plt.yticks(np.arange(0,len(ry)*nmax[1],len(ry))+(len(ry)/2),range(nmax[1]))
        plt.ylim(ymin=len(ry)*nmax[1],ymax=0)
        plt.title('Hermite Basis Functions (Cartesian)')
        plt.colorbar()

    if not (opts.savefig is None):
        plt.savefig(opts.savefig)
    else: plt.show()

