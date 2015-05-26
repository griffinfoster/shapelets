from pyrap.tables import table
import numpy as np
import os
import gc
import ephem
import re
import pyrap.quanta as qa
import pyrap.measures as pm
import rad2hmsdms


def Rotate(MS,ToRaDec):
    
    ra, dec = MS.radec
    ra1,dec1=ToRaDec

    # preparing rotation matrices (ORIGINAL PHASE DIRECTION)
    x = np.sin(ra)*np.cos(dec)
    y = np.cos(ra)*np.cos(dec)
    z = np.sin(dec)
    w = np.array([[x,y,z]]).T
    x = -np.sin(ra)*np.sin(dec)
    y = -np.cos(ra)*np.sin(dec)
    z = np.cos(dec)
    v = np.array([[x,y,z]]).T
    x = np.cos(ra)
    y = -np.sin(ra)
    z = 0
    u = np.array([[x,y,z]]).T
    T = np.concatenate([u,v,w], axis = -1 )


    TT=np.identity(3)
    x1 = np.sin(ra1)*np.cos(dec1)
    y1 = np.cos(ra1)*np.cos(dec1)
    z1 = np.sin(dec1)
    w1 = np.array([[x1,y1,z1]]).T
    x1 = -np.sin(ra1)*np.sin(dec1)
    y1 = -np.cos(ra1)*np.sin(dec1)
    z1 = np.cos(dec1)
    v1 = np.array([[x1,y1,z1]]).T
    x1 = np.cos(ra1)
    y1 = -np.sin(ra1)
    z1 = 0
    u1 = np.array([[x1,y1,z1]]).T

    Tshift = np.concatenate([u1,v1,w1], axis=-1)

    TT = np.dot(Tshift.T,T)

    uvw_all = MS.uvw
    uvw=uvw_all

    Phase=np.dot(np.dot((w-w1).T, T) , uvw.T)
    uvwNew=np.dot(uvw, TT.T)
    MS.uvw=uvwNew

    for chan in range(MS.NSPWChan):
        wavelength = MS.wavelength_chan.flatten()[chan]
        f = np.exp(Phase * 2 * np.pi * 1j/wavelength)
        for pol in range(4):
            MS.data[:,chan,pol]=MS.data[:,chan,pol] * f.reshape((MS.nrows,))


