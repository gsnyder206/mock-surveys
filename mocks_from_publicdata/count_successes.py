import glob
import os
import astropy.io.fits as fits
import numpy as np



fitslist=np.sort(np.asarray(glob.glob('*_F*.fits')))

for ff in fitslist:
    fd=fits.open(ff)[2].data
    print(ff,np.sum(fd['image_success']),fd.shape )
