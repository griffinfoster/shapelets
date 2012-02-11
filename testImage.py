#!/usr/bin/env python

import Image
import pylab as p

import fileio

fn='cyg_s15.png'
im=fileio.readImg(fn,gs=True)
print im.shape

p.imshow(im)
p.colorbar()
p.show()

#im=Image.open(fn)
#im0=im.convert("L")
#im0.show()

#p.subplot(221)
#p.imshow(im[:,:,0])
#p.subplot(222)
#p.imshow(im[:,:,1])
#p.subplot(223)
#p.imshow(im[:,:,2])
#p.subplot(224)
#p.imshow(im[:,:,3])
#p.show()

#im=Image.open(fn)

#im.show()

