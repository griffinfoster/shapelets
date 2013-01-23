shapelets
===

Created: 15.02.12 
Last Modified: 23.01.13  
Contact: griffin.foster@gmail.com  

Python module for fitting images (FITS,PNG,JPEG...) to shapelets, support for Cartesian and polar forms using
Hermite and Laguerre polynomials.  

Based on the shapelet framework developed in [Shapelets: I. A Method for Image Analysis](http://arxiv.org/abs/astro-ph/0105178) and the [IDL shapelet software](http://www.astro.caltech.edu/~rjm/shapelets/). 

Recommended Python Modules
===

* distutils 
* matplotlib/pylab 
* numpy 
* scipy 
* pyfits 
* PIL 
* cPickle (usually a standard install with python) 
* pyrap 

Install
===

sudo python setup.py install  

Usage
===

The scripts directory contains a number scripts for fitting.

* fitHermite.py: fits beta, and the centroid to use Hermite polynomials(Cartesian) 
* fitLaguerre.py: fits beta, and the centroid to use Laguerre polynomials(Polar) 
* fitBeta.py: fits only beta, uses a fixed centroid, both sets of polynomials can be fit for 
* solveShapelet.py: solve the shapelet decomposition for a given centroid, and beta (Polar or Cartesian) 
* plot*.py: scripts for plotting coefficient files, shapelet basis functions and images 

Examples
===

`scripts/fitBeta.py -r 170,335,195,309 --max data/lba_cyg_250.ms.CORRECTED_DATA.channel.1ch.fits`

