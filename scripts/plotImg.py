#!/usr/bin/env python
"""
Generic Image Plotter
"""

import sys,os
import pyfits as pf
import numpy as n
import pylab as p
import shapelets

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_IMAGE')
    o.set_description(__doc__)
    o.add_option('-r', '--region', dest='region', default=None,
        help='Region of image plot, (ymin,ymax,xmin,xmax), default: None')
    opts, args = o.parse_args(sys.argv[1:])

    fn=args[0]
    if fn.split('.')[-1].lower() == 'fits':
        im=shapelets.fileio.readFITS(fn)
    else:
        im=shapelets.fileio.readImg(fn)
    if not (opts.region is None):
        extent=map(int, opts.region.split(','))
        im=shapelets.img.selPxRange(im,extent)
    
    p.title('Image')
    p.imshow(im)
    p.colorbar()
    p.show()

