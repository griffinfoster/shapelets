
import numpy as np

def rad2hmsdms(rad,Type="dec",deg=False):
    
    if deg==False:
        deg=rad*180/np.pi
    else:
        deg=rad

    strsgn="+"
    if Type=="ra":
        deg/=15.
        strsgn=""

    if np.sign(deg)==-1.: strsgn="-"
    deg=np.abs(deg)
    degd=np.int(deg)
    degms=(deg-degd)*60.
    degm=np.int(degms)
    degs=((degms-degm)*60)
    return "%s%2.2i %2.2i %06.3f"%(strsgn,degd,degm,degs)