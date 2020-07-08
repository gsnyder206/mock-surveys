import cProfile
import pstats
import math
import string
import sys
import struct
import numpy as np
import cPickle
import asciitable
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import scipy.odr as odr
import glob
import os
import make_color_image
import make_fake_wht
import gzip
import tarfile
import shutil
import cosmocalc
import congrid
import astropy.io.ascii as ascii
import warnings
import subprocess
import photutils
import astropy
import astropy.cosmology
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.constants import G
from astropy.cosmology import WMAP7,z_at_value
from astropy.coordinates import SkyCoord
import copy
import medianstats_bootstrap as msbs


if __name__=="__main__":
    bb = pyfits.open('bbdata.fits')
    Npix = bb['STELLAR_MASS'].data.shape[0]
    if Npix != 4096:
        mtot = np.sum(bb['STELLAR_MASS'].data)
        mstar = congrid.congrid( bb['STELLAR_MASS'].data, (4096,4096),centre=True)
        mnew = np.sum(mstar)
        mstar = mstar*(mtot/mnew)
        print mnew, mtot

        zs = bb['WEIGHTED_Z'].data

        #sfrtot = np.sum(bb['SFR'].data)
        #sfr = congrid.congrid( bb['SFR'].data, (4096,4096),centre=True)
        #sfrnew = np.sum(sfr)
        #sfr = sfr*(sfrtot/sfrnew)

        z_f160 = congrid.congrid(zs[10,:,:], (4096,4096),centre=True)
        z_f444 = congrid.congrid(zs[18,:,:], (4096,4096),centre=True)
    else:
        mstar = bb['STELLAR_MASS'].data
        zs = bb['WEIGHTED_Z'].data
        #sfr = bb['SFR'].data
        z_f160 = zs[10,:,:]
        z_f444 = zs[18,:,:]

    ims = bb['PRIMARY'].data



    print np.log10(np.max(mstar))

    z1_hdu = pyfits.PrimaryHDU(z_f160)
    z1list = pyfits.HDUList([z1_hdu])
    z1list.writeto('REDSHIFTS_F160W.fits',clobber=True)

    z2_hdu = pyfits.PrimaryHDU(z_f444)
    z2list = pyfits.HDUList([z2_hdu])
    z2list.writeto('REDSHIFTS_F444W.fits',clobber=True)

    mstar_hdu = pyfits.PrimaryHDU(mstar)
    mstarlist = pyfits.HDUList([mstar_hdu])
    mstarlist.writeto('STELLARMASS.fits',clobber=True)

    #sfr_hdu = pyfits.PrimaryHDU(sfr)
    #sfrlist = pyfits.HDUList([sfr_hdu])
    #sfrlist.writeto('SFR.fits',clobber=True)

    if Npix != 4096:
        fils = bb['FILTERNAMES'].data.field('filter')
        for i,fil in enumerate(fils):
            filstr = fil.rstrip('.res')
            imoutf = filstr +'_regrid.fits'
            new_im = congrid.congrid(ims[i,:,:], (4096,4096), centre=True )
            im_hdu = pyfits.PrimaryHDU(new_im)
            imlist = pyfits.HDUList([im_hdu])
            imlist.writeto(imoutf,clobber=True)



    
    slices = ['0001','0002','0003','0004','0005','0006','0007','0008','0009']
    for sl in slices:
        filen = 'bbslice_'+sl+'.fits'
        if not os.path.lexists(filen):
            continue

        bb = pyfits.open(filen)

        mstar = bb['STELLAR_MASS'].data
        zs = bb['ZWEIGHTSLICE'].data
        ims = bb['SCI'].data

        z1_hdu = pyfits.PrimaryHDU(zs[10,:,:])
        z1list = pyfits.HDUList([z1_hdu])
        z1list.writeto('REDSHIFTS_F160W_'+sl+'.fits',clobber=True)

        z2_hdu = pyfits.PrimaryHDU(zs[18,:,:])
        z2list = pyfits.HDUList([z2_hdu])
        z2list.writeto('REDSHIFTS_F444W_'+sl+'.fits',clobber=True)

        mstar_hdu = pyfits.PrimaryHDU(mstar)
        mstarlist = pyfits.HDUList([mstar_hdu])
        mstarlist.writeto('STELLARMASS_'+sl+'.fits',clobber=True)



        fils = bb['FILTERS'].data.field('filter')
        

        
        for i,fil in enumerate(fils):
            
            filstr = fil.rstrip('.res')
            imoutf = filstr +'_'+ sl+'.fits'
            im_hdu = pyfits.PrimaryHDU(ims[i,:,:])
            imlist = pyfits.HDUList([im_hdu])
            imlist.writeto(imoutf,clobber=True)
        

        
