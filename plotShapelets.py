#!/usr/bin/env python
"""
Plot a set of shapelet basis functions
"""

import sys,os
import numpy as n
import pylab as p

#shapelet functions
import img, decomp, shapelet, fileio

#2x nmax
#2x beta

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
        help='Characteristic shapelet size, if using Hermite can use 2 values i.e. \'2.1,3.5\' default: 1.0')
    opts, args = o.parse_args(sys.argv[1:])

    nmax=opts.nmax.split(',')
    if len(nmax)==1: nmax=[int(nmax[0])+1,int(nmax[0])+1]
    else: nmax=[int(nmax[0])+1,int(nmax[1])+1]

    beta=opts.beta.split(',')
    if len(beta)==1: beta=[float(beta[0]),float(beta[0])]
    else: beta=[float(beta[0]),float(beta[1])]

    xlim=[-5,5]
    ylim=[-5,5]
    rx=n.arange(xlim[0],xlim[1],.1)
    ry=n.arange(ylim[0],ylim[1],.1)

    if opts.polar:
        nmax=nmax[0]
        beta=beta[0]
        print 'polar shapelets'
        fullImgReal=n.zeros((len(ry)*nmax*2,len(rx)*nmax))
        fullImgImag=n.zeros((len(ry)*nmax*2,len(rx)*nmax))
        yOffset=len(ry)*nmax
        r,th=shapelet.xy2rth(rx,ry)
        for nn in range(nmax):
            for mm in n.arange(-1*nn,nn+1):
                if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                    bf=shapelet.polarDimBasis(nn,mm,beta=beta)
                    bval=shapelet.computeBasisPolar(bf,r,th)
                    fullImgReal[mm*len(ry)+yOffset:(mm+1)*len(ry)+yOffset,nn*len(rx):(nn+1)*len(rx)]=bval.real
                    fullImgImag[mm*len(ry)+yOffset:(mm+1)*len(ry)+yOffset,nn*len(rx):(nn+1)*len(rx)]=bval.imag
        fig=p.figure()
        fig.subplots_adjust(wspace=0.3)
        
        p.subplot(121)
        p.imshow(fullImgReal)
        for nn in range(nmax):
            for mm in n.arange(-1*nn,nn+1):
                if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                    p.fill([nn*len(rx),(nn+1)*len(rx),(nn+1)*len(rx),nn*len(rx)],[mm*len(ry)+yOffset,mm*len(ry)+yOffset,(mm+1)*len(ry)+yOffset,(mm+1)*len(ry)+yOffset],fill=False)
        p.xlim(xmin=0,xmax=len(rx)*nmax)
        p.ylim(ymin=len(ry)*nmax*2,ymax=len(ry))
        p.xticks(n.arange(0,len(rx)*nmax,len(rx))+(len(rx)/2),range(nmax))
        p.yticks(n.arange(0,len(ry)*(nmax*2+1),len(ry)+1)+(len(ry)/2),nmax-n.arange(nmax*2+1))
        p.title('Real')
        p.colorbar()
        
        p.subplot(122)
        p.imshow(fullImgImag)
        for nn in range(nmax):
            for mm in n.arange(-1*nn,nn+1):
                if (nn%2==0 and mm%2==0) or (nn%2==1 and mm%2==1):
                    p.fill([nn*len(rx),(nn+1)*len(rx),(nn+1)*len(rx),nn*len(rx)],[mm*len(ry)+yOffset,mm*len(ry)+yOffset,(mm+1)*len(ry)+yOffset,(mm+1)*len(ry)+yOffset],fill=False)
        p.xlim(xmin=0,xmax=len(rx)*nmax)
        p.ylim(ymin=len(ry)*nmax*2,ymax=len(ry))
        p.xticks(n.arange(0,len(rx)*nmax,len(rx))+(len(rx)/2),range(nmax))
        p.yticks(n.arange(0,len(ry)*(nmax*2+1),len(ry)+1)+(len(ry)/2),nmax-n.arange(nmax*2+1))
        p.title('Imaginary')
        p.colorbar()
       
        p.suptitle('Lageurre Basis Functions (Polar)')
        p.show()
    
    else:
        print 'cartesian shapelets'
        fig=p.figure()
        fullImg=n.zeros((len(rx)*nmax[1],len(ry)*nmax[0]))
        for n0 in range(nmax[1]):
            for n1 in range(nmax[0]):
                bf=shapelet.dimBasis2d(n0,n1,beta=beta)
                bval=shapelet.computeBasis2d(bf,rx,ry)
                fullImg[n0*len(rx):(n0+1)*len(rx),n1*len(ry):(n1+1)*len(ry)]=bval
        
        p.imshow(fullImg)
        for n0 in range(nmax[1]):
            for n1 in range(nmax[0]):
                p.fill([n1*len(rx),(n1+1)*len(rx),(n1+1)*len(rx),n1*len(rx)],[n0*len(ry),n0*len(ry),(n0+1)*len(ry),(n0+1)*len(ry)],fill=False)
        p.xlim(xmin=0,xmax=len(rx)*nmax[0])
        p.xticks(n.arange(0,len(rx)*nmax[0],len(rx))+(len(rx)/2),range(nmax[0]))
        p.yticks(n.arange(0,len(ry)*nmax[1],len(ry))+(len(ry)/2),range(nmax[1]))
        p.ylim(ymin=len(ry)*nmax[1],ymax=0)
        p.title('Hermite Basis Functions (Cartesian)')
        p.colorbar()

        p.show()

