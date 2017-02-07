"""
A module interface to Measurement Set phase rotation functions
"""

import imp

try:
    imp.find_module('casacore')
    casacore = True
except ImportError:
    casacore = False

if casacore:
    import ClassMS, ModRotate # these modules require casacore

import ModColor, rad2hmsdms, reformat

