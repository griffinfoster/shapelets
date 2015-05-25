#!/usr/bin/env python
"""
Insert a shapelet coefficient file into the visibilities of a measurement set
"""

import sys
import numpy as np
import scipy.constants
from matplotlib import pyplot as plt
import shapelets
import pyrap.tables as pt
import distutils.dir_util

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] -c COEFF_FILE MS')
    o.set_description(__doc__)
    o.add_option('-c', '--coeff', dest='cfn', default=None,
        help='Shapelet coefficient file, default: None')
    o.add_option('-m', '--mode', dest='mode', default='add',
        help='\'replace\': Replace the visibilities, \'add\': add to the visibilites, \'subtract\': subtract from visibilites, default: add')
    o.add_option('-d', '--data_column', dest='data_column', default='DATA',
        help='Data column to take visibilities from/write to default: DATA')
    o.add_option('-o', '--outfile', dest='ofn', default=None,
        help='Output measurement set name, default: input file name + \'.shape\'')
    o.add_option('-X','--over', dest='overwrite', action="store_true", default=False,
        help='Over write original visibilities')
    opts, args = o.parse_args(sys.argv[1:])

    #load cofficients file
    if opts.cfn is None:
        sys.exit('error:missing shapelet coefficient file')
    d=shapelets.fileio.readCoeffs(opts.cfn)

    #scale beta to radians
    rotM=np.matrix([[np.cos(d['phi']),-1.*np.sin(d['phi'])],[np.sin(d['phi']),np.cos(d['phi'])]])
    rotDeltas=(np.pi/180.)*np.dot(rotM,np.array([d['dra'],d['ddec']])) #rotate delta RA and delta Dec, convert from degrees to radians
    betaRad=np.abs(np.multiply(np.array([d['beta']]),rotDeltas))

    #TODO: option: rephase to insert shapelet anywhere, default: at phase centre (?) is there a computationally easier way?
    #assume the shapelet is being insert at the phase centre (l,m)=(0,0)

    #generate basis functions
    print 'generating UV Basis Functions...',
    if d['mode'].startswith('herm'):
        #bfs=shapelets.decomp.genBasisFuncs([betaRad[0,0],betaRad[0,1]],d['norder'],d['phi'],fourier=True)
        print [betaRad[0,0],betaRad[0,1]]
        bfs=shapelets.decomp.genBasisFuncs([betaRad[0,0]*10.,betaRad[0,1]*10.],d['norder'],d['phi'],fourier=True) #TODO: my beta scaling looks to be off by a factor of 10!
        #bfs=shapelets.decomp.genBasisFuncs([.05,.05],d['norder'],d['phi'],fourier=True)
    elif d['mode'].startswith('lag'):
        bfs=shapelets.decomp.genPolarBasisFuncs([betaRad[0,0],betaRad[0,1]],d['norder'],d['phi'],fourier=True)
    print len(bfs), 'done'

    #load MS, get visibilites and u,v,w positions
    data_column=opts.data_column.upper()
    for fid,fn in enumerate(args):
        print 'working on: %s (%i of %i)'%(fn,fid+1,len(args))
        #ms=pt.table(fn,readonly=False)
        ms=pt.table(fn,readonly=True)
        uvw=ms.col('UVW').getcol() # [vis id, (u,v,w)]
        vis=ms.col(data_column).getcol() #[vis id, freq id, stokes id]
        #print uvw.shape, vis.shape
        #print uvw

        #gather channel frequency information
        sw=pt.table(fn+'/SPECTRAL_WINDOW')
        freqs=sw.col('CHAN_FREQ').getcol()
        cc=scipy.constants.c
        uu=uvw[:,0]
        vv=uvw[:,1]
        uu=np.reshape(uu,(uvw.shape[0],1))
        vv=np.reshape(vv,(uvw.shape[0],1))
        #convert u,v to units of wavelengths
        uu=(uu*freqs)/cc
        vv=(vv*freqs)/cc
        sw.close()

        #plt.plot(uu.flatten(),vv.flatten(),'r.')
        #plt.show()

        #evaulate basis functions at each u,v postion
        #TODO: with only real valued visibilities then there is mirror symetry in the image plane, where are the imaginary components?
        print 'Evaluating shapelet basis functions for all (u,v) positions'
        shapeVis=np.zeros_like(vis[:,:,0]).flatten() #only doing this for Stokes I
        for bfid,bf in enumerate(bfs):
            #print d['coeffs'][bfid]
            #print np.max(shapelets.shapelet.computeBasis2d(bf,uu.flatten(),vv.flatten()))
            shapeVis+=d['coeffs'][bfid]*shapelets.shapelet.computeBasis2d(bf,uu.flatten(),vv.flatten())
        #shapeVis+=1.*shapelets.shapelet.computeBasis2d(bfs[1],uu.flatten(),vv.flatten())
        shapeVis=np.reshape(shapeVis,uu.shape)
        print 'done'
        #print shapeVis

        #update visibilites
        #TODO: assuming fully linear polarization for the moment it is all going into the X receptor, this needs to be properly fixed
        print 'Updating visibilities (mode:%s)'%opts.mode
        if opts.mode.startswith('replace'):
            newVis=np.zeros_like(vis)
            newVis[:,:,0]=shapeVis
        elif opts.mode.startswith('add'):
            newVis=vis
            newVis[:,:,0]+=shapeVis
        elif opts.mode.startswith('sub'):
            newVis=vis
            newVis[:,:,0]-=shapeVis
        else:
            print 'Unknown MS update mode, the MS will be left unchanged'
        print 'done'

        ms.close()
        if opts.overwrite:
            oms=pt.table(fn,readonly=False)
            oms.putcol(data_column, newVis)
            oms.close()
        else:
            ofn=fn+'.shape'
            print '\twriting to:',ofn
            distutils.dir_util.copy_tree(fn,ofn)
            oms=pt.table(ofn,readonly=False)
            oms.putcol(data_column, newVis)
            oms.close()

