from setuptools import setup
import os, sys, glob

__version__ = '0.1'

setup(name = 'shapelets',
      version = __version__,
      description = 'Shapelet fitting and plotting',
      requires = ['pylab','numpy','pyfits','PIL','cPickle','scipy'],
      provides = ['shapelets'],
      package_dir = {'shapelets':'src'},
      packages = ['shapelets'],
      scripts=glob.glob('scripts/*.py'),
)
