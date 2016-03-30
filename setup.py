from setuptools import setup, find_packages
#from Cython.Build import cythonize
import numpy as np
import os, sys, glob

__version__ = '0.2'

setup(name = 'shapelets',
    version = __version__,
    description = 'Shapelet fitting and plotting',
    long_description = 'Shapelet fitting and plotting',
    author = 'Griffin Foster',
    author_email = 'griffin.foster@gmail.com',
    url = 'https://github.com/griffinfoster/shapelets',
    platforms = ['*nix'],
    license='GPL',
    requires = ['distutils', 'numpy', 'astropy', 'scipy', 'matplotlib', 'json'],
    provides = ['shapelets', 'shapelets.phs'],
    packages = ['shapelets', 'shapelets.phs'],
    #ext_modules = cythonize('shapelets/cshapelet.pyx', annotate=True),
    include_dirs = [np.get_include()],
    #scripts = glob.glob('scripts/*.py'),
    scripts = ['scripts/fitShapelet.py', 'scripts/insertShapelet.py', 'scripts/plotCoeffs.py', 'scripts/plotImg.py', 'scripts/plotShapelets.py', 'scripts/solveShapelet.py'],
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

