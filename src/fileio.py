"""
Read various image filetypes: FITS,png,jpeg
Write coeff file
"""

import pyfits as pf
import numpy as np
from PIL import Image
import cPickle as pickle


def readFITS(fn,hdr=False):
    """Read a FITS image file and returns a numpy array
    """
    hdulist=pf.open(fn)
    im=hdulist[0].data
    hdulist.close()
    if hdr:
        return im[0,0],getFITSInfo(fn)
    else: return im[0,0]

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

def getFITSInfo(fn):
    """Parse the FITS header for pointing and pixel size information
    return [RA,DEC], pixel resolution, pixel of [RA,DEC]
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
    ra=hdr['CRVAL1']
    dra=hdr['CDELT1']
    raPix=hdr['CRPIX1']
    dec=hdr['CRVAL2']
    ddec=hdr['CDELT2']
    decPix=hdr['CRPIX2']
    hdulist.close()
    return {'ra':ra,'dec':dec,'dra':dra,'ddec':ddec,'raPix':raPix,'decPix':decPix}

def writeHermiteCoeffs(fn,coeffs,xc,size,beta,norder,pos=[0.,0.,0.,0.],mode='hermite',info=''):
    """Write hermite coeffs and meta data to a pickle file
    fn: output file name
    coeffs: set of coefficients for Hermite polynomials
    xc: center position
    size: size in pixels of image
    beta: characteristic beta values
    norder: order of polynomials
    mode: basis function mode
    pos: 4 element array of RA,DEC and sale size from FITS header
    info: extra metadata space
    """
    d={ 'coeffs':coeffs,
        'mode':mode,
        'xc':xc,
        'size':size,
        'beta':beta,
        'norder':norder,
        'ra':pos[0],
        'dec':pos[1],
        'dra':pos[2],
        'ddec':pos[2],
        'info': info }
    fh=open(fn,'wb')
    pickle.dump(d,fh)
    fh.close()

def readHermiteCoeffs(fn):
    """Read a binary coeff file and return a dictionary of values
    fn: filename to read
    """
    fh=open(fn,'rb')
    d=pickle.load(fh)
    fh.close()
    return d

def writeLageurreCoeffs(fn,coeffs,xc,size,beta,norder,pos=[0.,0.,0.,0.],mode='laguerre',info=''):
    """Write Lageurre coeffs and meta data to a pickle file
    fn: output file name
    coeffs: set of coefficients for Lageurre polynomials
    xc: center position
    size: size in pixels of image
    beta: characteristic beta value
    norder: max order of polynomials
    mode: basis function mode
    pos: 4 element array of RA,DEC and sale size from FITS header
    info: extra metadata space
    """
    d={ 'coeffs':coeffs,
        'mode':mode,
        'xc':xc,
        'size':size,
        'beta':beta,
        'norder':norder,
        'ra':pos[0],
        'dec':pos[1],
        'dra':pos[2],
        'ddec':pos[2],
        'info': info }
    fh=open(fn,'wb')
    pickle.dump(d,fh)
    fh.close()

def readLageurreCoeffs(fn):
    """Read a binary coeff file and return a dictionary of values
    fn: filename to read
    """
    fh=open(fn,'rb')
    d=pickle.load(fh)
    fh.close()
    return d

def readCoeffs(fn):
    """At present readLageurreCoeffs and readHermiteCoeffs do the same operations"""
    return readHermiteCoeffs(fn)

if __name__ == "__main__":
    print "fileio"
