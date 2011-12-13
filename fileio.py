"""
Read various filetypes: FITS

ToDo:
    read png/jpeg images
"""

import pyfits as pf
import pylab as p

import img

if __name__ == "__main__":
    fn='../lba_cyg_20101105_200.ms.CORRECTED_DATA.channel.1ch.fits'
    hdulist=pf.open(fn)
    #print hdulist.info()
    #print hdulist[0].header
    im=hdulist[0].data
    #im=img.selPxRange(im[0,0],(100,400,150,420))
    im=im[0,0]
    #print img.flux(im)
    #centroid=img.centroid(im)
    centroid=img.centroid(im,region=[200,300,200,300])
    #print centroid
    #noisemap=img.estimateNoiseMap(im,masks=[[136,176,82,128]],sigma=3.,maxiter=10)
    noisemap=img.estimateNoiseMap(im,sigma=3.,maxiter=10)
    p.imshow(im)
    #p.imshow(noisemap)
    p.text(centroid[0],centroid[1], '+')
    p.colorbar()
    p.show()
