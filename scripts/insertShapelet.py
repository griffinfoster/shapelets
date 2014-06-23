#!/usr/bin/env python
"""
Replace measurement set vis data with a shapelet model
"""

import sys
import distutils.dir_util
import numpy as np
import shapelets
import pyrap.tables as pt

import pylab as p

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
    print 'RA: %f\t DEC: %f\t DRA: %f\t DDEC: %f\t'%(d['ra'],d['dec'],d['dra'],d['ddec'])
    beta=[d['beta'][0]*d['dra'],d['beta'][1]*d['ddec']]
    print 'BETA: (%f,%f)\t INVERSE: (%f,%f)'%(beta[0],beta[1],1./beta[0],1./beta[1])
    print 'generating UV Basis Functions...',
    sys.stdout.flush()
    if d['mode'].startswith('herm'):
        bfs=shapelets.shapelet.ftHermiteBasis([(np.pi/180.)*beta[0],(np.pi/180.)*beta[1]],d['norder'])
    elif d['mode'].startswith('lag'):
        bfs=shapelets.shapelet.ftLaguerreBasis(d['dra'],d['norder'])
    print 'done'

    #image->uv transform limits
    dra=d['dra']*np.pi/180.
    ddec=d['ddec']*np.pi/180.
    rx=np.array(range(0,d['size'][0]),dtype=float)-d['xc'][0]
    ry=np.array(range(0,d['size'][1]),dtype=float)-d['xc'][1]
    usize=(float(len(ry))/float(len(rx)))*1./(np.abs(rx[0]-rx[1])*dra) #!!!
    vsize=(float(len(rx))/float(len(ry)))*1./(np.abs(ry[0]-ry[1])*ddec) #!!!
    #print usize,vsize
    ures=(float(len(ry))/float(len(rx)))*1./(np.abs(rx[0]-rx[-1])*dra) #!!!
    vres=(float(len(rx))/float(len(ry)))*1./(np.abs(ry[0]-ry[-1])*ddec) #!!!
    #print ures,vres
    umin=ures*(rx[0]+.5)
    vmin=vres*(ry[0]+.5)
    uu=np.arange(0,usize+ures,ures)+umin
    vv=np.arange(0,vsize+vres,vres)+vmin
    uRange=[uu[0],uu[-1]]
    vRange=[vv[0],vv[-1]]

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
        u=np.reshape(u,(uvwData.shape[0],1))
        v=np.reshape(v,(uvwData.shape[0],1))
        #convert u,v to units of wavelengths
        u=(u*freqs.T)/cc
        v=(v*freqs.T)/cc
        sw.close()

        #bvals=[]
        #for bf in bfs:
        #    bvals.append(shapelets.shapelet.computeBasis2d(bf,uu,vv).flatten())
        #bm=np.array(bvals)
        #mdl_ft=shapelets.img.constructModel(bm.transpose(),d['coeffs'],[len(uu),len(vv)])
        #p.imshow(np.abs(mdl_ft))
        #p.show()
        
        #filter out UV samples
        uvFilter=((u>uRange[0]) & (u<uRange[1]) & ((v>vRange[0]) & (v<vRange[1])))
        uvIndex=np.argwhere(uvFilter)
        u0=u[uvFilter]
        v0=v[uvFilter]

        #compute the UV complex correlation from coefficients and basis functions
        print '\tcomputing new UV visibilities...',
        sys.stdout.flush()
        if d['mode'].startswith('herm'):
            sVis=shapelets.uv.computeHermiteUV(bfs,d['coeffs'],u0,v0)
        elif d['mode'].startswith('lag'):
            sVis=shapelets.uv.computeLaguerreUV(bfs,d['coeffs'],u,v)
        print 'done'

        #uv coverage plot
        p.scatter(u0,v0,c=np.abs(sVis)/np.max(np.abs(sVis)),edgecolor='none')
        p.show()

        vis=ms.col(data_column)
        visData=vis.getcol()
        visShape=visData.shape
        newVis=np.zeros_like(visData)
        xxVis=newVis[:,:,0]
        xxVis[uvIndex[:,0],0]=sVis
        newVis[:,:,0]=xxVis     #only writing to the XX data
        
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

