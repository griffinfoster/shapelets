from setuptools import setup, find_packages
#from Cython.Build import cythonize
import numpy as np
import os, sys, glob

__version__ = '0.2'

setup(name = 'shapelets',
    version = __version__,
    description = 'Shapelet fitting and plotting',
    long_description = 'Shapelet fitting and plotting',
    author='Griffin Foster',
    author_email='griffin.foster@gmail.com',
    url='',
    platforms=['*nix'],
    license='GPL',
    #requires = ['distutils','numpy>=1.8.2','pyfits>=3.2','scipy>=0.13.3','matplotlib','json','pywcs'],
    requires = ['distutils','numpy','pyfits','scipy','matplotlib','json','pywcs'],
    dependency_links = ['http://www.stsci.edu/resources/software_hardware/pyfits'],
    provides = ['shapelets', 'shapelets.phs'],
    packages = ['shapelets', 'shapelets.phs'],
    #ext_modules = cythonize('shapelets/cshapelet.pyx', annotate=True),
    include_dirs=[np.get_include()],
    #scripts=glob.glob('scripts/*.py'),
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
)

