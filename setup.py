from setuptools import setup
from Cython.Build import cythonize
import numpy as np
import cython_gsl
import os, sys, glob

__version__ = '0.1'

setup(name = 'shapelets',
      version = __version__,
      description = 'Shapelet fitting and plotting',
      long_description = 'Shapelet fitting and plotting',
      author='Griffin Foster',
      author_email='griffin.foster@gmail.com',
      url='',
      platforms=['*nix'],
      license='GPL',
      requires = ['pyrap','distutils','pylab','numpy','pyfits','PIL','cPickle','scipy'],
      provides = ['shapelets'],
      package_dir = {'shapelets':'src'},
      packages = ['shapelets'],
      ext_modules = cythonize('shapelets/cshapelet.pyx', annotate=True),
      include_dirs=[np.get_include(), cython_gsl.get_include()],
      libraries=[('gsl',{'sources':['shapelets/cshapelet.pyx']}),
                 ('gslcblas', {'sources': ['shapelets/cshapelet.pyx']})],#cython_gsl.get_libraries(),
      scripts=glob.glob('scripts/*.py'),
)
