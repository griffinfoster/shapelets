from setuptools import setup
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
      scripts=glob.glob('scripts/*.py'),
)
