import astropy
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import numpy as np
import matplotlib
import matplotlib.pyplot as pyplot
import make_color_image
import astropy.units as u
from astropy.cosmology import WMAP7,z_at_value
import glob
import sys
import os

def pd_rgb(angle='0',alph=3.0e-5,Q=5.0):

    r=0.75*1.0e-35*((np.load('image.1.0.angle'+angle+'.npz'))['image'])
    g=1.0*1.0e-35*((np.load('image.0.5.angle'+angle+'.npz'))['image'])
    b=2.5*1.0e-35*((np.load('image.0.3.angle'+angle+'.npz'))['image'])

    savefile='powderday_rgb_angle'+angle+'.pdf'
    print(savefile)

    f=pyplot.figure(figsize=(4.0,4.0),dpi=600)
    pyplot.subplots_adjust(wspace=0.0,hspace=0.0,top=1.0,right=1.0,left=0.00,bottom=0.00)

    ax=f.add_subplot(1,1,1)
    ax.set_xticks([]) ; ax.set_yticks([])

    rgbthing = make_color_image.make_interactive_nasa(b,g,r,alph,Q)

    ax.imshow(rgbthing,interpolation='nearest',aspect='auto',origin='lower')

    f.savefig(savefile,dpi=600)
    pyplot.close(f)

    return
