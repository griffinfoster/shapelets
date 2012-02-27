#!/usr/bin/env python
"""
Replace measurement set vis data with a shapelet model
"""

import sys,os
import distutils.dir_util
import numpy as n
import shapelets

import pyrap.tables as pt

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] -c COEFF_FILE MS')
    o.set_description(__doc__)
    o.add_option('-c', '--coeff', dest='cfn', default=None,
        help='Shapelet coefficient file, default: None')
    o.add_option('-m', '--mode', dest='mode', default='replace',
        help='\'replace\': Replace the visibilities, \'add\': add to the visibilites, default: replace')
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
    
    #generate inverse basis functions
    print 'RA: %f\t DEC: %f\t BETA: (%f,%f)\t INVERSE: (%f,%f)'%(d['ra'],d['dec'],d['dra'],d['ddec'],1./d['dra'],1./d['ddec'])
    print 'generating UV correlations...',
    sys.stdout.flush()
    if d['mode'].startswith('herm'):
        bfs=shapelets.shapelet.ftHermiteBasis([(n.pi/180.)*d['dra'],(n.pi/180.)*d['ddec']],d['norder'])
        #bfs=shapelets.shapelet.ftHermiteBasis([d['dra'],d['ddec']],d['norder'])
    elif d['mode'].startswith('lag'):
        bfs=shapelets.shapelet.ftLaguerreBasis(d['dra'],d['norder'])
    print 'done'
    
    data_column=opts.data_column.upper()
    for fn in args:
        print 'working on:',fn
        ms=pt.table(fn,readonly=False)
        uvw=ms.col('UVW')
        uvwData=uvw.getcol()
        
        #gather channel frequency information
        sw=pt.table(fn+'/SPECTRAL_WINDOW')
        chan_freqs=sw.col('CHAN_FREQ')
        freqs=chan_freqs.getcol()
        cc=299792458.0
        u=uvwData[:,0]
        v=uvwData[:,1]
        u=n.reshape(u,(uvwData.shape[0],1))
        v=n.reshape(v,(uvwData.shape[0],1))
        u=(u*freqs.T)/cc
        v=(v*freqs.T)/cc
        sw.close()

        #compute the UV complex correlation from coefficients and basis functions
        print '\tcomputing new UV visibilities...',
        sys.stdout.flush()
        if d['mode'].startswith('herm'):
            sVis=shapelets.uv.computeHermiteUV(bfs,d['coeffs'],u,v)
        elif d['mode'].startswith('lag'):
            sVis=shapelets.uv.computeLaguerreUV(bfs,d['coeffs'],u,v)
        print 'done'
        
        vis=ms.col(data_column)
        visData=vis.getcol()
        visShape=visData.shape
        newVis=n.zeros_like(visData)
        newVis[:,:,0]=sVis  #only writing to the XX data
        
        if opts.mode.startswith('add'): newVis+=visData
        if opts.overwrite:
            ms.putcol(data_column, newVis)
            ms.close()
        else:
            ms.close()
            ofn=fn+'.shape'
            print '\twriting to:',ofn
            distutils.dir_util.copy_tree(fn,ofn)
            ms=pt.table(ofn,readonly=False)
            vis=ms.col(data_column)
            ms.putcol(data_column, newVis)

