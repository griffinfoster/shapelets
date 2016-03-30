shapelets
===

Contact: griffin.foster@gmail.com  

A python module for fitting and decomposing images (FITS, PNG, JPEG...) into shapelet coefficients, support for Cartesian and polar forms using Hermite and Laguerre polynomials.  

Based on the shapelet framework developed in [Shapelets: I. A Method for Image Analysis](http://arxiv.org/abs/astro-ph/0105178) and the [IDL shapelet software](http://www.astro.caltech.edu/~rjm/shapelets/). 

#### Required Python Modules
 
* matplotlib
* numpy
* scipy
* astropy
* json

#### Optional Python Modules

* [python-casacore](https://github.com/casacore/python-casacore)
* PIL 

#### Install

While developing it is useful to do a developer install:

```
sudo python setup.py develop
```

Otherwise, the standard install will work:

```
sudo python setup.py install  
```

#### Usage

The scripts directory contains a number scripts for plotting, decomposing, and fitting.

* plotShapelets.py : plot a grid of shapelet basis functions
* plotImg.py : load a FITS or image file, useful for determing coordinates to apply shapelet decomposition
* plotCoeffs.py : given a shapelet coefficient file, plot the modelled source and coefficients
* solveShapelet.py : given a set of parameters and an image, decompose the image into shapelet coefficients
* fitShapelet.py : fit parameters to minimize the chi^2 difference between shapelet model and image
* insertShapelet.py : insert a shapelet coefficient set into a Measurement Set (requires python-casacore)

#### Examples

```
plotShapelets.py -n 4 -p -b 1.0,1.0,0.44

plotImg.py -r 1010,1117,947,1030 ../data/N6251_test.fits
plotImg.py ../data/zen.2455819.69771.uvcRREM_briggs-dirty.fits

solveShapelet.py -r 1010,1117,947,1030 -N 891,1257,600,840 ../data/N6251_test.fits -n 15 -x 49,52 --beta=6.,2.5 --phi=-0.448243 -m cart
solveShapelet.py -r 489,530,489,527 -N 436,561,405,487 ../data/zen.2455819.69771.uvcRREM_briggs-dirty.fits -n 10 -m cart

fitShapelet.py -r 1010,1117,947,1030 -N 891,1257,600,840 ../../data/N6251_test.fits -n 8 -x 49,52 --init_beta=6.,2.5 --init_phi=-0.448243 -m cart
fitShapelet.py -r 489,530,489,527 -N 436,561,405,487 ../../data/zen.2455819.69771.uvcRREM_briggs-dirty.fits -n 8 -x 20.,22. --init_beta=2.788068,2.336974 --init_phi=-1.046530 -m cart -B 10

plotCoeffs.py ../data/solve/N6251_cart.pkl

```
