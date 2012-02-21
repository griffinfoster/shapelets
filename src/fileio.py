"""
Read various image filetypes: FITS,png,jpeg
Write coeff file
"""

import pyfits as pf
import numpy as n
from PIL import Image
import cPickle as pickle

import img

def readFITS(fn):
    """Read a FITS image file and returns a numpy array
    """
    hdulist=pf.open(fn)
    im=hdulist[0].data
    return im[0,0]

def readImg(fn,gs=False):
    """Read an image file using PIL libraries and return a numpy array
    valid types(tested other may work): png, jpeg
    gs: return a greyscale image
    """
    im=Image.open(fn)
    if gs: im=im.convert("L")
    return n.asarray(im)

def writeHermiteCoeffs(fn,coeffs,xc,size,beta,norder,mode='hermite',info=''):
    """Write hermite coeffs and meta data to a pickle file
    fn: output file name
    coeffs: set of coefficients for Hermite polynomials
    xc: center position
    size: size in pixels of image
    beta: characteristic beta values
    norder: order of polynomials
    mode: basis function mode
    info: extra metadata space
    """
    d={ 'coeffs':coeffs,
        'mode':mode,
        'xc':xc,
        'size':size,
        'beta':beta,
        'norder':norder,
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

def writeLageurreCoeffs(fn,coeffs,xc,size,beta,norder,mode='laguerre',info=''):
    """Write Lageurre coeffs and meta data to a pickle file
    fn: output file name
    coeffs: set of coefficients for Lageurre polynomials
    xc: center position
    size: size in pixels of image
    beta: characteristic beta value
    norder: max order of polynomials
    mode: basis function mode
    info: extra metadata space
    """
    d={ 'coeffs':coeffs,
        'mode':mode,
        'xc':xc,
        'size':size,
        'beta':beta,
        'norder':norder,
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
