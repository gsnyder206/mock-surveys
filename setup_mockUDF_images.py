#import ez_galaxy
import math
import string
import sys
import struct
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.colors as pycolors
#import utils
import numpy as np
#import array
#import astLib.astStats as astStats
import cPickle
import asciitable
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import pyfits
import glob
import os
import make_color_image
import make_fake_wht
import gzip
import tarfile
import shutil
#import psfmatch
from pygoods import *
from sextutils import *
import cosmocalc
#from parse_hydroART_morphs import *
import numpy.ma as ma
import congrid


sq_arcsec_per_sr = 42545170296.0
c = 3.0e8

def do_psf(initial_image, PSF_use, image_pixsize, psf_pixsize):

    orig_size = initial_image.shape[0]

    size = orig_size*image_pixsize/psf_pixsize

    print orig_size, size

    binned_image = congrid.congrid(initial_image, (size,size),method='linear')

    print binned_image.shape

    conv_im = scipy.signal.convolve2d(binned_image,PSF_use/np.sum(PSF_use),mode='same',boundary='wrap')

    print conv_im.shape

    out_im = congrid.congrid(conv_im,(orig_size,orig_size),method='linear')

    print out_im.shape

    return out_im





if __name__=="__main__":

    #initial_image_file = '/Users/gsnyder/Documents/Projects/mockUDFs/FITS_FILES/fromOdyssey/YXZ/DECEMBER30_1e7_new_ism_both/bbslice_0006.fits'
    #initial_image_file = '/Users/gsnyder/Documents/Projects/mockUDFs/FITS_FILES/fromOdyssey/YXZ/DECEMBER30_1e7_new_ism_both/bbdata.fits'

    #initial_image_dir = '/Users/gsnyder/Documents/Projects/mockUDFs/FITS_FILES/fromOdyssey/ZYX/DECEMBER30_1e7_new_ism_both'  #FIELDC
    #initial_image_dir = '/Users/gsnyder/Documents/Projects/mockUDFs/FITS_FILES/fromOdyssey/YXZ/DECEMBER30_1e7_new_ism_both' #FIELDA
    initial_image_dir = '/Users/gsnyder/Documents/Projects/mockUDFs/FITS_FILES/fromOdyssey/XYZ/DECEMBER30_1e7_new_ism_both'  #FIELDB

    bbs1 = os.path.join(initial_image_dir,'bbslice_0001.fits')
    bbs2 = os.path.join(initial_image_dir,'bbslice_0002.fits')
    bbs3 = os.path.join(initial_image_dir,'bbslice_0003.fits')
    bbs4 = os.path.join(initial_image_dir,'bbslice_0004.fits')
    bbs5 = os.path.join(initial_image_dir,'bbslice_0005.fits')
    bbs6 = os.path.join(initial_image_dir,'bbslice_0006.fits')
    bbs7 = os.path.join(initial_image_dir,'bbslice_0007.fits')
    bbs8 = os.path.join(initial_image_dir,'bbslice_0008.fits')
    bbs9 = os.path.join(initial_image_dir,'bbslice_0009.fits')

    initial_image = pyfits.open(bbs1)[0].data +pyfits.open(bbs2)[0].data +pyfits.open(bbs3)[0].data + pyfits.open(bbs4)[0].data + pyfits.open(bbs5)[0].data + pyfits.open(bbs6)[0].data + pyfits.open(bbs7)[0].data + pyfits.open(bbs8)[0].data + pyfits.open(bbs9)[0].data

    #initial_image = pyfits.open('/Users/gsnyder/Documents/Projects/mockUDFs/FITS_FILES/fromOdyssey/ZYX/DECEMBER30_1e7_new_ism_both/broadband_0065.fits')[12].data

    print pyfits.open(bbs1)[0].header.cards, pyfits.open(bbs9)[0].header.cards
    print initial_image.shape

    use_psf=True
    label = 'FIELDB'


    #initial_image_J = (pyfits.open(initial_image_file)[0].data)[8]
    #initial_image_z = (pyfits.open(initial_image_file)[0].data)[6]

    #efflam = 1.60
    #to_nJy_per_Sr = (1.0e9)*(1.0e14)*(efflam**2)/c   #((pixscale/206265.0)^2)*
    #initial_image_nJySr = to_nJy_per_Sr*initial_image
    #to_nJy_per_Sr = (1.0e9)*(1.0e14)*(1.25**2)/c   #((pixscale/206265.0)^2)*
    #initial_image_J_nJySr = to_nJy_per_Sr*initial_image_J    
    #to_nJy_per_Sr = (1.0e9)*(1.0e14)*(0.90**2)/c   #((pixscale/206265.0)^2)*
    #initial_image_z_nJySr = to_nJy_per_Sr*initial_image_z

    #f160, F125, F850 == index 10, 8, 6
    #slice 0006 === z=2-3

    ko=40
    kf=60
    #PSF_h = (pyfits.open('/Users/gsnyder/Documents/Projects/HydroART_Morphology/NONPARAM_CATALOGS/HST_PSFS/F160W_rebin.fits')[0].data)[ko:kf,ko:kf]
    wfc3_scale = 0.13
    acs_scale = 0.05
    ##????

    jpsf_dir = os.path.expandvars('$HOME/Dropbox/Professional/Illustris_Morphology/AndersonPSFs')

    fPSF_h = os.path.join(jpsf_dir,'PSFSTD_WFC3IR_F160W.fits')
    fPSF_140 = os.path.join(jpsf_dir,'PSFSTD_WFC3IR_F140W.fits')
    fPSF_j = os.path.join(jpsf_dir,'PSFSTD_WFC3IR_F125W.fits')
    fPSF_105 = os.path.join(jpsf_dir,'PSFSTD_WFC3IR_F105W.fits')
    fPSF_336 = os.path.join(jpsf_dir,'PSFSTD_WFC3UV_F336W.fits')

    fPSF_850 = os.path.join(jpsf_dir,'PSFSTD_ACSWFC_F850L_SM3.fits')
    fPSF_775 = os.path.join(jpsf_dir,'PSFSTD_ACSWFC_F775W_SM3.fits')
    fPSF_814 = os.path.join(jpsf_dir,'PSFSTD_ACSWFC_F814W.fits')
    fPSF_606 = os.path.join(jpsf_dir,'PSFSTD_ACSWFC_F606W.fits')
    fPSF_435 = os.path.join(jpsf_dir,'PSFSTD_ACSWFC_F435W.fits')

    wfc3_psf_pix = wfc3_scale/4.0
    acs_psf_pix = acs_scale/4.0
    #nircam??

    PSF_h = ((pyfits.open(fPSF_h)[0]).data)[0,ko:kf,ko:kf]
    PSF_140 = ((pyfits.open(fPSF_140)[0]).data)[0,ko:kf,ko:kf]
    PSF_j = ((pyfits.open(fPSF_j)[0]).data)[0,ko:kf,ko:kf]
    PSF_105 = ((pyfits.open(fPSF_105)[0]).data)[0,ko:kf,ko:kf]
    #PSF_336 = ((pyfits.open(fPSF_336)[0]).data)[0,ko:kf,ko:kf]

    PSF_850 = ((pyfits.open(fPSF_850)[0]).data)[0,ko:kf,ko:kf]
    PSF_775 = ((pyfits.open(fPSF_775)[0]).data)[0,ko:kf,ko:kf]
    PSF_814 = ((pyfits.open(fPSF_814)[0]).data)[0,ko:kf,ko:kf]
    PSF_606 = ((pyfits.open(fPSF_606)[0]).data)[0,ko:kf,ko:kf]
    PSF_435 = ((pyfits.open(fPSF_435)[0]).data)[0,ko:kf,ko:kf]



    pixel_arcsec = 0.04161785759037315  #arcsec
    pixel_Sr = (pixel_arcsec**2)/sq_arcsec_per_sr
    FWHM_arcsec = 0.115
    Npix = (initial_image.shape)[-1]

    #exit()


    output_dir_base = '/Users/gsnyder/Documents/Projects/mockUDFs/SE_RUNS/'+label
    if not os.path.lexists(output_dir_base):
        os.mkdir(output_dir_base)


    cex_directory = "/Users/gsnyder/Documents/PythonCode/cex/"
    cex_files = np.asarray(glob.glob(os.path.join(cex_directory,"*")))

    maglims = [27.0] #mag/arcsec^2
    maglabs = ['SB27']

    whtdir = ''
    filters = ['F160W','F125W','F850LP']#,'F775W','F606W','F435W','F814W','F105W','F140W']
    indices = [10,8,6,4,3,2,5,7,9]
    efflams = [1.60,1.25,0.90,0.775,0.606,0.435,0.814,1.05,1.40]
    PSFlist = [PSF_h,PSF_j,PSF_850,PSF_775,PSF_606,PSF_435,PSF_814,PSF_105,PSF_140]
    PSF_pixlist = [wfc3_psf_pix,wfc3_psf_pix,acs_psf_pix,acs_psf_pix,acs_psf_pix,acs_psf_pix,acs_psf_pix,wfc3_psf_pix,wfc3_psf_pix]

    i=0
    for maglim in maglims:
        depthdir = os.path.join(output_dir_base,label+'_'+maglabs[i])
        if not os.path.lexists(depthdir):
            os.mkdir(depthdir)

        sigma_nJy = (2.0**(-0.5))*((1.0e9)*(3631.0/5.0)*10.0**(-0.4*maglim))*pixel_arcsec*(3.0*FWHM_arcsec)

        for cexf in cex_files:
            shutil.copy(cexf,depthdir)


        j=0
        for fil in filters:
            imbase = label+'_'+maglabs[i]+'_'+fil

            imname = os.path.join(depthdir,imbase+'.fits')

            to_nJy_per_Sr = (1.0e9)*(1.0e14)*(efflams[j]**2)/c   #((pixscale/206265.0)^2)*
            to_microJy_per_arcsec = (1.0e-3)*(to_nJy_per_Sr/sq_arcsec_per_sr)

            #initial_image_nJySr = to_nJy_per_Sr*initial_image[indices[j]]
            #AB_zeropoint = -2.5*np.log10(pixel_Sr) - 2.5*(-9.0) + 2.5*np.log10(3631.0)
            AB_zeropoint = -2.5*np.log10(pixel_arcsec**2) - 2.5*(-6.0) + 2.5*np.log10(3631.0)
            initial_image_uJyAs = to_microJy_per_arcsec*initial_image[indices[j]]

            if use_psf==True:

                PSF_use = PSFlist[j]
                PSF_pix = PSF_pixlist[j]
                print PSF_use.shape, PSF_pix, j
                #convolved_image_nJySr = do_psf(initial_image_nJySr,PSF_use,pixel_arcsec,PSF_pix)
                convolved_image_uJyAs = do_psf(initial_image_uJyAs,PSF_use,pixel_arcsec,PSF_pix)
                #print convolved_image_nJySr.shape
                #j=j+1 ; continue
            else:
                #convolved_image_nJySr = initial_image_nJySr
                convolved_image_uJyAs = initial_image_uJyAs



            #final_image = (convolved_image_nJySr + (sigma_nJy/pixel_Sr)*np.random.randn(Npix,Npix))
            final_image = (convolved_image_uJyAs + ((sigma_nJy*1.0e-3)/(pixel_arcsec**2))*np.random.randn(Npix,Npix))


            HDU = pyfits.PrimaryHDU(final_image)
            (HDU.header)['BUNIT'] = 'uJy/Arcsec' #'nJy/Sr'
            (HDU.header)['PIXSCALE'] = (pixel_arcsec, 'arcsec')
            (HDU.header)['ABZP'] = AB_zeropoint
            (HDU.header)['RMS'] = (sigma_nJy*1.0e-3)/(pixel_arcsec**2) #(sigma_nJy/pixel_Sr)
            (HDU.header)['RMSnJy'] = sigma_nJy

            HDUlist = pyfits.HDUList([HDU])
            HDUlist.writeto(imname,clobber=True) 

            wht_name_link = os.path.join(depthdir, imbase+'_wht.fits')
            flag_name_link = os.path.join(depthdir, imbase+'_flag.fits')
            invvar_name_link = os.path.join(depthdir, imbase+'_invvar.fits')

            whtcheck = os.path.join(whtdir, 'wht_file_4096')
            flagcheck = os.path.join(whtdir, 'flag_file_4096')

            valcheck = os.path.join(whtdir, 'val_file_4096_'+fil+'_'+maglabs[i])

            if not os.path.lexists(valcheck):
                whtfile, flagfile, varfile = make_fake_wht.make_maps(Nacross = int(Npix), value=(1.0/(sigma_nJy/pixel_Sr))**2, name = '_'+fil+'_'+maglabs[i])  #value=?
            else:
                whtfile=whtcheck
                flagfile=flagcheck
                varfile=valcheck


            if not os.path.lexists(wht_name_link):
                os.symlink(whtfile, wht_name_link)
            if not os.path.lexists(flag_name_link):
                os.symlink(flagfile, flag_name_link)
            if not os.path.lexists(invvar_name_link):
                os.symlink(varfile, invvar_name_link)


  

            j = j+1




        setupfile = os.path.join(depthdir,'SETUP_sex_mockUDF.txt')
        setupobj = open(setupfile,'r+')
        newlines = []
        for line in setupobj:
            swapsb = line.replace('SB29',maglabs[i])
            swaplabel = swapsb.replace('FIELDA',label)
            newlines.append(swaplabel)
        setupobj.close()
        newobj = open(setupfile,'w')
        for line in newlines:
            newobj.write(line)
        newobj.close()
                  
        i = i+1
