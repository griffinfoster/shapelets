#!/usr/bin/env python
"""
Insert a shapelet coefficient file into the visibilities of a measurement set
"""

import sys
import numpy as np
import scipy.constants
#from matplotlib import pyplot as plt
import shapelets
import shapelets.phs
import casacore.tables as tbls
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
        help='A flux scaling factor. This is useful if the model was derived from an averaged wideband image and the coefficients need to be rescaled to reverse the averaging, in this situation the rescaling value should be the number of channels imaged. A value of 1 will cause no rescaling Default: Number of channels in the MS')

    #rephasing options
    o.add_option('-p','--phs',dest='phsMode', default='centre',
        help='Placement of source by phase rotating the visibilities, options are centre: the model is inserted at the phase centre no phase rotation is performed ; source: visibilities are rotated to make the centroid of the model the phase centre, then rotated back to the original phase centre ; manual: a manually input (RA,dec) is rotated to and the source inserted. default: centre')
    o.add_option('--ra', dest='ra', default=None,
        help='RA (in radians unless flag set) to phase data to')
    o.add_option('--dec', dest='dec', default=None,
        help='DEC (in radians unless flag set) to phase data to')
    o.add_option('--deg',dest='deg_flag',action='store_true',
        help='Use degrees instead of radians')
    o.add_option('--str',dest='str_flag', action='store_true',
        help='Use hh:mm:ss.sss format for RA and dd:mm:ss.sss format for DEC')

    o.add_option('-b','--beta',dest='beta',default=None,
        help='Override beta (in pixels), can be a 1, 2, or 3 comma seperated values')
    opts, args = o.parse_args(sys.argv[1:])


    #load coefficient file
    if opts.cfn is None:
        sys.exit('error:missing shapelet coefficient file')
    d=shapelets.fileio.readCoeffs(opts.cfn)

    #override beta if input as option
    if opts.beta is None:
        beta=d['beta']
        phi=d['phi']
    else:
        betaList=map(float,opts.beta.split(','))
        if len(betaList)==1:
            beta=[betaList[0],betaList[0]]
            phi=d['phi']
        elif len(betaList)==2:
            beta=betaList
            phi=d['phi']
        elif len(betaList)==3:
            beta=[betaList[0],betaList[1]]
            phi=betaList[2]

    ##scale beta to radians
    #rotM=np.matrix([[np.cos(d['phi']),-1.*np.sin(d['phi'])],[np.sin(d['phi']),np.cos(d['phi'])]])
    #temp0=np.dot(rotM,np.array([(np.pi/180.)*d['dra'],0.]))
    #temp1=np.dot(rotM,np.array([0.,(np.pi/180.)*d['ddec']]))
    #rotDeltas=np.array([np.sqrt(temp0[0,0]**2.+temp0[0,1]**2.), np.sqrt(temp1[0,0]**2.+temp1[0,1]**2.)])
    #betaRad=(2.*np.pi)*rotDeltas*np.array(d['beta']) #there is some mathematical convention issue, something about wavenumber and wavelength...should sort out a clear answer
    #phi=d['phi']

    #scale beta to radians
    rotM=np.matrix([[np.cos(phi),-1.*np.sin(phi)],[np.sin(phi),np.cos(phi)]])
    temp0=np.dot(rotM,np.array([(np.pi/180.)*d['dra'],0.]))
    temp1=np.dot(rotM,np.array([0.,(np.pi/180.)*d['ddec']]))
    rotDeltas=np.array([np.sqrt(temp0[0,0]**2.+temp0[0,1]**2.), np.sqrt(temp1[0,0]**2.+temp1[0,1]**2.)])
    betaRad=(2.*np.pi)*rotDeltas*np.array(beta) #there is some mathematical convention issue, something about wavenumber and wavelength...should sort out a clear answer

    #generate basis functions
    print 'generating UV Basis Functions...',
    if d['mode'].startswith('herm'):
        bfs=shapelets.decomp.genBasisFuncs([betaRad[0],betaRad[1]],d['norder'],phi,fourier=True)
    elif d['mode'].startswith('lag'):
        bfs=shapelets.decomp.genPolarBasisFuncs([betaRad[0],betaRad[1]],d['norder'],phi,fourier=True)
    print len(bfs), 'done'

    #parse where to insert the shapelet model
    if opts.phsMode.startswith('manual'):
        if opts.ra is None or opts.dec is None:
            print 'ERROR: RA or DEC not set'
            exit(1)

        if opts.deg_flag:
            ra=np.pi*float(opts.ra)/180.
            dec=np.pi*float(opts.dec)/180.
        elif opts.str_flag:
            raArr=map(float,opts.ra.split(':'))
            ra=15.*(raArr[0]+raArr[1]/60.+raArr[2]/3600.)
            ra=np.pi*ra/180.
            decArr=map(float,opts.dec.split(':'))
            dec=np.abs(decArr[0])+decArr[1]/60.+decArr[2]/3600.
            if decArr[0] < 0.: dec*=-1.
            dec=np.pi*dec/180.
        else:
            ra=float(opts.ra)
            dec=float(opts.dec)
        phsRotate=True
    elif opts.phsMode.startswith('source'):
        ra=np.pi*d['ra']/180.
        dec=np.pi*d['dec']/180.
        phsRotate=True
    else:
        phsRotate=False

    if phsRotate:
        print 'Inserting shapelet model at RA=%f Dec=%f (mode=%s)'%(ra,dec,opts.phsMode)
    else: 
        print 'Inserting shapelet model at phase centre'

    #load MS, get visibilites and u,v,w positions
    data_column=opts.data_column.upper()
    for fid,fn in enumerate(args):

        if opts.overwrite:
            ofn=fn
        else:
            ofn=fn+'.shape'
            distutils.dir_util.copy_tree(fn,ofn)
        print 'working on: %s (%i of %i)'%(ofn,fid+1,len(args))

        #get the initial phase centre
        field=pt.table(ofn+'/FIELD')
        phaseCentre=field.col('PHASE_DIR').getcol()[0,0]
        field.close()
        if phsRotate:
            print 'Phasing to RA: %1.10f, DEC: %1.10f'%(ra,dec)
            MS=shapelets.phs.ClassMS.ClassMS(ofn,Col=data_column)
            MS.RotateMS((ra,dec))
            MS.SaveVis(Col=data_column)
            print 'done'

        MS=tbls.table(ofn,readonly=True)
        uvw=MS.col('UVW').getcol() # [vis id, (u,v,w)]
        vis=MS.col(data_column).getcol() #[vis id, freq id, stokes id]
        MS.close()

        #gather channel frequency information
        SW=tbls.table(ofn+'/SPECTRAL_WINDOW')
        freqs=SW.col('CHAN_FREQ').getcol()
        cc=scipy.constants.c
        uu=uvw[:,0]
        vv=uvw[:,1]
        uu=np.reshape(uu,(uvw.shape[0],1))
        vv=np.reshape(vv,(uvw.shape[0],1))
        #convert u,v to units of wavelengths
        uu=(uu*freqs)/cc
        vv=(vv*freqs)/cc
        SW.close()

        #rescaling value
        #flux scaling is too low if the model was derived form a wideband averaged image, in that case the correct rescaling is the number of channels
        if opts.rescale is None:
            rescale=float(freqs.shape[1])
        else: rescale=float(opts.rescale)
        print 'Rescale factor: %f'%rescale

        #evaulate basis functions at each u,v postion
        print 'Evaluating shapelet basis functions for all (u,v) positions'
        shapeVis=np.zeros_like(vis[:,:,0]).flatten() #only doing this for Stokes I
        for bfid,bf in enumerate(bfs):
            #TODO: visibilites seem to be rotated by 90 degrees
            #TODO: i think this line is correct, but the rotation and mirroring leads me to use the line below -> shapeVis+=d['coeffs'][bfid]*shapelets.shapelet.computeBasis2d(bf,uu.flatten(),vv.flatten())
            #TODO: this could be at least partially due to the sign of the FITS ddec and dra
            shapeVis+=d['coeffs'][bfid]*shapelets.shapelet.computeBasis2d(bf,-1.*vv.flatten(),-1.*uu.flatten())
        #shapeVis+=1.*shapelets.shapelet.computeBasis2d(bfs[1],uu.flatten(),vv.flatten())
        shapeVis=rescale*np.reshape(shapeVis,uu.shape)
        print 'done'

        #update visibilites
        #TODO: assuming unpolarized source for the moment, this needs to be fixed properly
        print 'Updating visibilities (mode:%s)'%opts.mode
        if opts.mode.startswith('replace'):
            newVis=np.zeros_like(vis)
            newVis[:,:,0]=shapeVis/2. #XX
            newVis[:,:,3]=shapeVis/2. #YY
        elif opts.mode.startswith('add'):
            newVis=vis
            newVis[:,:,0]+=shapeVis/2. #XX
            newVis[:,:,3]+=shapeVis/2. #YY
        elif opts.mode.startswith('sub'):
            newVis=vis
            newVis[:,:,0]-=shapeVis/2. #XX
            newVis[:,:,3]-=shapeVis/2. #YY
        else:
            print 'Unknown MS update mode, the MS will be left unchanged'
        print 'done'
        
        print 'writing to:',ofn
        MS=tbls.table(ofn,readonly=False)
        MS.putcol(data_column, newVis)
        MS.close()
        print 'done'

        if phsRotate:
            print 'Phase rotating back to original phase centre...'
            print 'Phasing to RA: %1.10f, DEC: %1.10f'%(phaseCentre[0],phaseCentre[1])
            MS=shapelets.phs.ClassMS.ClassMS(ofn,Col=data_column)
            MS.RotateMS((phaseCentre[0],phaseCentre[1]))
            MS.SaveVis(Col=data_column)
            print 'done'

