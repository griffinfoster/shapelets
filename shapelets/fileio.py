"""
Read various image filetypes: FITS,png,jpeg
Write coeff file
"""

import pyfits as pf
import numpy as np
import cPickle as pickle
import sys
import json
import pywcs

#optional packages:
try:
    from PIL import Image
except ImportError: pass

def readFITS(fn,hdr=False):
    """Read a FITS image file and returns a numpy array (only Stokes I or the first Stokes index)
    """
    hdulist=pf.open(fn)
    im=hdulist[0].data
    #image data format: [frequency, polarization, dec, ra]
    hdulist.close()
    h=getFITSInfo(fn)
    im=im[0,0]
    #orient the image to the same way it would look in a normal FITS viewer
    if h['dra'] > 0: im=np.fliplr(im)
    if h['ddec'] > 0: im=np.flipud(im)
    if hdr: return im,h
    else: return im

def getFITSInfo(fn):
    """Parse the FITS header for pointing and pixel size information
    return [RA,DEC], pixel resolution, pixel of [RA,DEC]
    generates a WCS instance for converting between sky and pixels
    """
    hdulist=pf.open(fn)
    hdr=hdulist[0].header
    #CTYPE1: RA---[PROJ], projection SIN/TAN/ARC
    #CRVAL1: reference RA position in degrees
    #CRPIX1: location of reference pixel
    #CDELT1: delta RA/pixel size in degrees
    #CTYPE2: DEC--[PROJ], projection SIN/TAN/ARC
    #CRVAL2: reference DEC position in degrees
    #CRPIX2: location of reference pixel
    #CDELT2: delta DEC/pixel size in degrees
    #BMAJ: major axis of PSF FWHM ellipse (degrees)
    #BMIN: minor axis of PSF FWHM ellipse (degrees)
    #BPA: rotation angle of PSF FWHM ellipse (degrees)
    ra=hdr['CRVAL1']
    dra=hdr['CDELT1']
    raPix=hdr['CRPIX1']
    dec=hdr['CRVAL2']
    ddec=hdr['CDELT2']
    decPix=hdr['CRPIX2']
    bmaj=hdr['BMAJ']
    bmin=hdr['BMIN']
    bpa=hdr['BPA']
    #Generate a WCS structure, using the normal method creates errors due to header formating
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crval = [ra,dec]
    wcs.wcs.crpix = [raPix,decPix]
    wcs.wcs.cdelt = [dra,ddec]
    wcs.wcs.ctype = ["RA---SIN", "DEC--SIN"]

    hdulist.close()
    return {'ra':ra,'dec':dec,'dra':dra,'ddec':ddec,'raPix':raPix,'decPix':decPix, 'wcs':wcs, 'bmaj':bmaj, 'bmin':bmin, 'bpa':bpa}

def readImg(fn,gs=False):
    """Read an image file using PIL libraries and return a numpy array
    valid types(tested other may work): png, jpeg
    gs: return a greyscale image
    """
    im=Image.open(fn)
    if gs: im=im.convert("L")
    return np.asarray(im)

def readArrayPkl(fn):
    """Read a pickle file, expected format is a NxM numpy array
    """
    fh=open(fn,'rb')
    im=pickle.load(fh)
    fh.close()
    return im

def writeHermiteCoeffs(fn,coeffs,xc,size,beta,phi,norder,pos=[0.,0.,0.,0.],mode='hermite',info='',fmt='pkl'):
    """Write hermite coeffs and meta data to a pickle file
    fn: output file name
    coeffs: set of coefficients for Hermite polynomials
    xc: center position
    size: size in pixels of image
    beta: characteristic beta values
    phi: rotation angle
    norder: order of polynomials
    mode: basis function mode
    pos: 4 element array of RA,DEC and scale size from FITS header (degrees)
    info: extra metadata space
    fmt: output formats supported: pkl (pickle), json
    """
    if type(norder) is int: norder=[norder,norder]
    d={ 'coeffs':coeffs,
        'mode':'hermite',
        'xc':xc,
        'size':size,
        'beta':beta,
        'phi':phi,
        'norder':norder,
        'ra':pos[0],
        'dec':pos[1],
        'dra':pos[2],
        'ddec':pos[3],
        'info': info }
    if fmt.startswith('pkl'):
        fh=open(fn,'wb')
        pickle.dump(d,fh)
        fh.close()
    elif fmt.startswith('json'):
        d['coeffs']=coeffs.tolist()
        fh=open(fn,'w')
        json.dump(d,fh,sort_keys=True)
        fh.close()
    else:
        print 'Unknown format, no file written'

def readHermiteCoeffs(fn):
    """Read a binary coeff file and return a dictionary of values
    fn: filename to read
    """
    if fn.endswith('.pkl'):
        fh=open(fn,'rb')
        d=pickle.load(fh)
        fh.close()
        return d
    elif fn.endswith('.json'):
        fh=open(fn,'r')
        d=json.load(fh)
        coeffs=np.array(d.pop('coeffs'))
        d['coeffs']=coeffs
        fh.close()
        return d
    else:
        print 'Unknown file type'
        return np.nan

def writeLageurreCoeffs(fn,coeffs,xc,size,beta,phi,norder,pos=[0.,0.,0.,0.],mode='laguerre',info='',fmt='pkl'):
    """Write Lageurre coeffs and meta data to a pickle file
    fn: output file name
    coeffs: set of coefficients for Lageurre polynomials
    xc: center position
    size: size in pixels of image
    beta: characteristic beta values
    phi: rotation angle
    norder: max order of polynomials
    mode: basis function mode
    pos: 4 element array of RA,DEC and scale size from FITS header (degrees)
    info: extra metadata space
    fmt: output formats supported: pkl (pickle), json
    """
    if type(norder) is int: norder=[norder,norder]
    d={ 'mode':'lageurre',
        'xc':xc,
        'size':size,
        'beta':beta,
        'phi':phi,
        'norder':norder,
        'ra':pos[0],
        'dec':pos[1],
        'dra':pos[2],
        'ddec':pos[3],
        'info': info }
    if fmt.startswith('pkl'):
        d['coeffs']=coeffs
        fh=open(fn,'wb')
        pickle.dump(d,fh)
        fh.close()
    elif fmt.startswith('json'):
        d['rcoeffs']=coeffs.real.tolist()
        d['icoeffs']=coeffs.imag.tolist()
        fh=open(fn,'w')
        json.dump(d,fh,sort_keys=True)
        fh.close()
    else:
        print 'Unknown format, no file written'

def readLageurreCoeffs(fn):
    """Read a binary coeff file and return a dictionary of values
    fn: filename to read
    """
    if fn.endswith('.pkl'):
        fh=open(fn,'rb')
        d=pickle.load(fh)
        fh.close()
        return d
    elif fn.endswith('.json'):
        fh=open(fn,'r')
        d=json.load(fh)
        rcoeffs=np.array(d.pop('rcoeffs'))
        icoeffs=np.array(d.pop('icoeffs'))
        d['coeffs']=rcoeffs + 1j*icoeffs
        fh.close()
        return d
    else:
        print 'Unknown file type'
        return np.nan

#TODO: a general coeff file reader
def readCoeffs(fn):
    """Determine a coefficient file based on file extension, call the correct reader
    """
    if fn.endswith('.pkl'):
        fh=open(fn,'rb')
        d=pickle.load(fh)
        fh.close()
        return d
    elif fn.endswith('.json'):
        fh=open(fn,'r')
        d=json.load(fh)
        fh.close()
        if d['mode'].startswith('hermite'): return readHermiteCoeffs(fn)
        elif d['mode'].startswith('lageurre'): return readLageurreCoeffs(fn)


if __name__ == "__main__":

    print '============================================'
    print 'Testing fileio module:'
    print '============================================'
    tc=0
    te=0
    
    #read in a FITS file
    tc+=1
    try:
        im,hdr=readFITS('../data/N6251_test.fits',hdr=True)
        print 'FITS header:', hdr
        print 'FITS image shape:', im.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #read in a PNG file
    tc+=1
    try:
        im=readImg('../data/N6251_test.png',gs=True)
        print 'PNG shape:', im.shape
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    #load pre computed shapelet coeffs (pkl)
    #save pre computed shapelet coeffs (pkl)
    #load pre computed shapelet coeffs (json)
    #save pre computed shapelet coeffs (json)
    tc+=1
    try:
        shapeDict=readHermiteCoeffs('../data/testHermite.pkl')
        print shapeDict.keys()
        
        writeHermiteCoeffs('testWriteHermite.pkl',shapeDict['coeffs'],shapeDict['xc'],shapeDict['size'],shapeDict['beta'],shapeDict['phi'],shapeDict['norder'],pos=[shapeDict['ra'],shapeDict['dec'],shapeDict['dra'],shapeDict['ddec']],info=shapeDict['info'],fmt='pkl')
        shapeDict0=readHermiteCoeffs('testWriteHermite.pkl')
        print shapeDict0.keys()

        writeHermiteCoeffs('testWriteHermite.json',shapeDict['coeffs'],shapeDict['xc'],shapeDict['size'],shapeDict['beta'],shapeDict['phi'],shapeDict['norder'],pos=[shapeDict['ra'],shapeDict['dec'],shapeDict['dra'],shapeDict['ddec']],info=shapeDict['info'],fmt='json')
        shapeDict0=readHermiteCoeffs('testWriteHermite.json')
        print shapeDict0.keys()

        shapeDict=readLageurreCoeffs('../data/testLageurre.pkl')
        print shapeDict.keys()
        
        writeLageurreCoeffs('testWriteLageurre.pkl',shapeDict['coeffs'],shapeDict['xc'],shapeDict['size'],shapeDict['beta'],shapeDict['phi'],shapeDict['norder'],pos=[shapeDict['ra'],shapeDict['dec'],shapeDict['dra'],shapeDict['ddec']],info=shapeDict['info'],fmt='pkl')
        shapeDict0=readLageurreCoeffs('testWriteLageurre.pkl')
        print shapeDict0.keys()
        
        writeLageurreCoeffs('testWriteLageurre.json',shapeDict['coeffs'],shapeDict['xc'],shapeDict['size'],shapeDict['beta'],shapeDict['phi'],shapeDict['norder'],pos=[shapeDict['ra'],shapeDict['dec'],shapeDict['dra'],shapeDict['ddec']],info=shapeDict['info'],fmt='json')
        shapeDict0=readLageurreCoeffs('testWriteLageurre.json')
        print shapeDict0.keys()
    except:
        print 'Test failed (%i):'%tc, sys.exc_info()[0]
        te+=1

    print '============================================'
    print '%i of %i tests succeeded'%(tc-te,tc)
    print '============================================'

