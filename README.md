shapelets
===

Created: 15.02.12  
Last Modified: 22.05.15  
Contact: griffin.foster@gmail.com  

A python module for fitting and decomposing images (FITS,PNG,JPEG...) into shapelet coefficients, support for Cartesian and polar forms using Hermite and Laguerre polynomials.  

Based on the shapelet framework developed in [Shapelets: I. A Method for Image Analysis](http://arxiv.org/abs/astro-ph/0105178) and the [IDL shapelet software](http://www.astro.caltech.edu/~rjm/shapelets/). 

#### Required Python Modules
===

* distutils (usually a standard install with python)
* cPickle (usually a standard install with python) 
* matplotlib 
* numpy 
* scipy 
* pyfits 
* json

#### Recommended Python Modules
===

* PIL 
* pyrap 

#### Install
===

While developing it is useful to do a developer install:

```
sudo python setup.py develop
```

Otherwise, the standard install will work:

```
sudo python setup.py install  
```

#### Usage
===

The scripts directory contains a number scripts for plotting, decomposing, and fitting.

* plotShapelets.py : plot a grid of shapelet basis functions
* plotImg.py : load a FITS or image file, useful for determing coordinates to apply shapelet decomposition
* plotCoeffs.py : given a shapelet coefficient file, plot the modelled source and coefficients
* solveShapelet.py : given a set of parameters and an image, decompose the image into shapelet coefficients
* fitShapelet.py : fit parameters to minimize the chi^2 difference between shapelet model and image

#### Examples
===

```
plotShapelets.py -n 4 -p -b 1.0,1.0,0.44

solveShapelet.py ../data/N6251_test.fits -r 1028,1097,1025,1074 -N 967,1067,972,1026 -n 20 -m polar
solveShapelet.py ../data/N6251_test.fits -r 1028,1097,1025,1074 -N 967,1067,972,1026 -n 20 -b 3.,3. -p 0.448266

fitShapelet.py ../data/N6251_test.fits -r 1028,1097,1025,1074 -N 967,1067,972,1026 -n 8 -m polar -p 0.486636 -b 4.,4. --set_phi --set_xc
fitShapelet.py ../data/N6251_test.fits -r 1028,1097,1025,1074 -N 967,1067,972,1026 -n 8 -p 0.486636 -b 4.,4. --set_phi --set_xc

```

