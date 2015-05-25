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
    o.add_option('-S','--scale',dest='rescale', default=None,
        help='A flux scaling factor. This is useful if the model was derived from an averaged wideband image and the coefficients need to be rescaled to reverse the averaging. A value of 1 will cause no rescaling Default: Number of channels in the MS')
    opts, args = o.parse_args(sys.argv[1:])

    #load coefficient file
    if opts.cfn is None:
        sys.exit('error:missing shapelet coefficient file')
    d=shapelets.fileio.readCoeffs(opts.cfn)

    #scale beta to radians
    rotM=np.matrix([[np.cos(d['phi']),-1.*np.sin(d['phi'])],[np.sin(d['phi']),np.cos(d['phi'])]])
    rotDeltas=(np.pi/180.)*np.dot(rotM,np.array([d['dra'],d['ddec']])) #rotate delta RA and delta Dec, convert from degrees to radians
    betaRad=(2.*np.pi)*np.abs(np.multiply(np.array([d['beta']]),rotDeltas)) #there is some mathematical convention issue, something about wavenumber and wavelength...should sort out a clear answer

    print d['beta'],(np.pi/180.)*np.array([d['dra'],d['ddec']]),rotDeltas,betaRad

    #TODO: option: rephase to insert shapelet anywhere, default: at phase centre (?) is there a computationally easier way?
    #assume the shapelet is being insert at the phase centre (l,m)=(0,0)

    #generate basis functions
    print 'generating UV Basis Functions...',
    if d['mode'].startswith('herm'):
        bfs=shapelets.decomp.genBasisFuncs([betaRad[0,0],betaRad[0,1]],d['norder'],d['phi'],fourier=True)
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

        #rescaling value
        #flux scaling is too low if the model was derived form a wideband averaged image, in that case the correct rescaling is the number of channels
        if opts.rescale is None:
            rescale=float(freqs.shape[1])
        else: rescale=float(opts.rescale)
        print 'Rescale factor: %f'%rescale

        #plt.plot(uu.flatten(),vv.flatten(),'r.')
        #plt.show()

        #evaulate basis functions at each u,v postion
        print 'Evaluating shapelet basis functions for all (u,v) positions'
        shapeVis=np.zeros_like(vis[:,:,0]).flatten() #only doing this for Stokes I
        for bfid,bf in enumerate(bfs):
            #TODO: visibilites seem to be rotated by 90 degrees
            #TODO: i think this line is correct, but the rotation and mirroring leads me to use the line below -> shapeVis+=d['coeffs'][bfid]*shapelets.shapelet.computeBasis2d(bf,uu.flatten(),vv.flatten())
            shapeVis+=rescale*d['coeffs'][bfid]*shapelets.shapelet.computeBasis2d(bf,vv.flatten(),-1.*uu.flatten())
        #shapeVis+=1.*shapelets.shapelet.computeBasis2d(bfs[1],uu.flatten(),vv.flatten())
        shapeVis=np.reshape(shapeVis,uu.shape)
        print 'done'
        #print shapeVis

        #update visibilites
        #TODO: assuming fully linear polarization for the moment, so all the flux is going into the X receptor, this needs to be properly fixed
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

