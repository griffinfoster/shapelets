from setuptools import setup, find_packages
#from Cython.Build import cythonize
import numpy as np
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
      requires = ['pyrap','distutils','pylab','numpy','pyfits','PIL','cPickle','scipy','json'],
      provides = ['shapelets'],
      packages = ['shapelets'],
      #package_dir = {'shapelets':'src'},
      #packages = find_packages(),
      #ext_modules = cythonize('shapelets/cshapelet.pyx', annotate=True),
      include_dirs=[np.get_include()],
      #scripts=glob.glob('scripts/*.py'),
)
