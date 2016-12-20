
####
#### Name:    hudfimagescript.py
#### Author:  Greg Snyder   gsnyder@stsci.edu
#### Purpose:  The final step.  Converts output of add_all_broadband_flux.py into PDF image files.
#### Disclaimer:  This code is provided AS-IS with absolutely NO warranty.
####           It is largely meant as a guide rather than ideal code.  
####           It can and should be replaced with other lightcone generation techniques. 
####           I make no claims about the immediate usability of this code.
####           That said, I am happy to field questions and discuss issues 
####           related to this code. And also to help use it.
#### License:  ?
#### Credit:   Users should cite Snyder et al. (2014): "Mock UDFs from Hydro Sims with Sunrise"
####           AND ALSO   Kitzbichler & White 2007, Henriques et al. 2013, Overzier et al. 2013, etc.
####           These papers were used as source material to create the default lightcone algorithm.
####  


import math
import string
import sys
import struct
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import asciitable
import scipy.ndimage
import scipy.ndimage.interpolation
import scipy.stats as ss
import scipy as sp
import pyfits
import glob
import os
import make_color_image




def realdatascript():


    folder = '/n/scratch2/hernquist_lab/gsnyder/XDF'
    outdir = '/n/home09/gsnyder/cosmo_stuff/XDF'

    f435w = pyfits.open(os.path.join(folder,'hlsp_xdf_hst_acswfc-60mas_hudf_f435w_v1_sci.fits'))

    print f435w.info()

    print f435w[0].header.ascardlist()


    f606w = pyfits.open(os.path.join(folder,'hlsp_xdf_hst_acswfc-60mas_hudf_f606w_v1_sci.fits'))
    f775w = pyfits.open(os.path.join(folder,'hlsp_xdf_hst_acswfc-60mas_hudf_f775w_v1_sci.fits'))
    f814w = pyfits.open(os.path.join(folder,'hlsp_xdf_hst_acswfc-60mas_hudf_f814w_v1_sci.fits'))


    f850lp = pyfits.open(os.path.join(folder, 'hlsp_xdf_hst_acswfc-60mas_hudf_f850lp_v1_sci.fits'))
    f160w = pyfits.open(os.path.join(folder, 'hlsp_xdf_hst_wfc3ir-60mas_hudf_f160w_v1_sci.fits'))


    f105w = pyfits.open(os.path.join(folder, 'hlsp_xdf_hst_wfc3ir-60mas_hudf_f105w_v1_sci.fits'))
    f125w = pyfits.open(os.path.join(folder, 'hlsp_xdf_hst_wfc3ir-60mas_hudf_f125w_v1_sci.fits'))
    f140w = pyfits.open(os.path.join(folder, 'hlsp_xdf_hst_wfc3ir-60mas_hudf_f140w_v1_sci.fits'))


    print f160w.info()
    print f160w[0].header.ascardlist()

    pixelscale = 0.06 # arcseconds.   0.06 = 60 mas



    #AB zeropoints from XDF website http://archive.stsci.edu/prepds/xdf/#dataproducts

    zp_f435w = 25.68
    zp_f606w = 26.51
    zp_f775w = 25.69
    zp_f814w = 25.94
    zp_f850lp = 24.87

    zp_f105w = 26.27
    zp_f125w = 26.23
    zp_f140w = 26.45
    zp_f160w = 25.94

    #most likely:  AB = -2.5 * log10{ IMAGE[electrons/s] } + ZP
    #  f/3631Jy = 10.0^{0.4*2.5log[IM]+ZP} = IM*10.0^ZP

    #also some PHOTFLAMs from a table that differs by ~2-3% from those above

    pf_f435w = 3.184049e-19  # erg/cm^2/Ang/count
    pf_f606w = 7.862481e-20
    pf_f775w = 9.969696e-20
    pf_f814w = 7.033186e-20
    pf_f850lp = 1.542299e-19

    #WFC3 appear to be a "perfect" match, though website says 2%
    #http://www.stsci.edu/hst/wfc3/phot_zp_lbn

    pf_f105w = 3.0386e-20
    pf_f125w = 2.2484e-20
    pf_f140w = 1.4737e-20
    pf_f160w = 1.9276e-20

    #all assume a count rate in e/s

    outputblock = np.ndarray(shape=(8,5250,5250))

    XDF_F435W = (4.25e17 * f435w[0].data * pf_f435w * ( pixelscale**(-2.0) ) )
    XDF_F606W = (4.25e17 * f606w[0].data * pf_f606w * ( pixelscale**(-2.0) ) )
    XDF_F775W = (4.25e17 * f775w[0].data * pf_f775w * ( pixelscale**(-2.0) ) )
    XDF_F814W = (4.25e17 * f814w[0].data * pf_f814w * ( pixelscale**(-2.0) ) )
    XDF_F850LP = (4.25e17 * f850lp[0].data * pf_f850lp * ( pixelscale**(-2.0) ) )
    
    XDF_F105W = (4.25e17 * f105w[0].data * pf_f105w * ( pixelscale**(-2.0) ) )
    XDF_F125W = (4.25e17 * f125w[0].data * pf_f125w * ( pixelscale**(-2.0) ) )
    XDF_F140W = (4.25e17 * f140w[0].data * pf_f140w * ( pixelscale**(-2.0) ) )
    XDF_F160W = (4.25e17 * f160w[0].data * pf_f160w * ( pixelscale**(-2.0) ) )



    F435W_conversion = 4.25e17 * pf_f435w * ( pixelscale**(-2.0) )
    F850LP_conversion = 4.25e17 * pf_f850lp * ( pixelscale**(-2.0) )
    F160W_conversion = 4.25e17 * pf_f160w * ( pixelscale**(-2.0) )

    print F435W_conversion, F850LP_conversion, F160W_conversion
    print F435W_conversion*0.0013, F850LP_conversion*0.0005, F160W_conversion*0.0009

    outputfilename = 'HST_XDF2012_SunriseUnits_60mas.fits'
    primhdu = pyfits.PrimaryHDU(data=None) ; primhdu.header.update('FILETYPE', 'HSTXDFIMAGES')


    f435w_hdu = pyfits.ImageHDU(XDF_F435W)
    f435w_hdu.header.update('EXTNAME','F435W')

    f606w_hdu = pyfits.ImageHDU(XDF_F606W)
    f606w_hdu.header.update('EXTNAME','F606W')

    f775w_hdu = pyfits.ImageHDU(XDF_F775W)
    f775w_hdu.header.update('EXTNAME','F775W')

    f814w_hdu = pyfits.ImageHDU(XDF_F814W)
    f814w_hdu.header.update('EXTNAME','F814W')

    f850lp_hdu = pyfits.ImageHDU(XDF_F850LP)
    f850lp_hdu.header.update('EXTNAME','F850LP')

    f105w_hdu = pyfits.ImageHDU(XDF_F105W)
    f105w_hdu.header.update('EXTNAME','F105W')
    f125w_hdu = pyfits.ImageHDU(XDF_F125W)
    f125w_hdu.header.update('EXTNAME','F125W')
    f140w_hdu = pyfits.ImageHDU(XDF_F140W)
    f140w_hdu.header.update('EXTNAME','F140W')
    f160w_hdu = pyfits.ImageHDU(XDF_F160W)
    f160w_hdu.header.update('EXTNAME','F160W')


    #unitlist = pyfits.HDUList([primhdu, f435w_hdu,f606w_hdu,f775w_hdu,f814w_hdu,f850lp_hdu,f105w_hdu,f125w_hdu,f140w_hdu,f160w_hdu])
    #unitlist.writeto(outputfilename, clobber=True)


    #From notebook analysis Fall 2012
    #Flux = 4.25e17 * (IM [electrons/s] ) * (pf_ [photflam]) * (pixelscale/arcsec)^-2    [W/m^2/m/Sr]




    XDF_F435W = (4.25e17 * f435w[0].data * pf_f435w * ( pixelscale**(-2.0) ) )[1204:4045,1204:4045]
    XDF_F850LP = (4.25e17 * f850lp[0].data * pf_f850lp * ( pixelscale**(-2.0) ) )[1204:4045,1204:4045] 
    XDF_F160W = (4.25e17 * f160w[0].data * pf_f160w * ( pixelscale**(-2.0) ) )[1204:4045,1204:4045]


    print (f435w[0].data).shape
    Rotated_f435w = scipy.ndimage.interpolation.rotate(f435w[0].data,42)
    Rotated_f850lp = scipy.ndimage.interpolation.rotate(f850lp[0].data,42)
    Rotated_f160w = scipy.ndimage.interpolation.rotate(f160w[0].data,42)

    print Rotated_f435w.shape
    XDF_F435W_rotated = (4.25e17 * Rotated_f435w * pf_f435w * ( pixelscale**(-2.0) ) )[2200:5041,2200:5041]   #4045-1204 = 2841
    XDF_F850LP_rotated = (4.25e17 * Rotated_f850lp * pf_f850lp * ( pixelscale**(-2.0) ) )[2200:5041,2200:5041]   #4045-1204 = 2841
    XDF_F160W_rotated = (4.25e17 * Rotated_f160w * pf_f160w * ( pixelscale**(-2.0) ) )[2200:5041,2200:5041]   #4045-1204 = 2841
    print XDF_F850LP_rotated.shape, XDF_F850LP.shape

    #result = make_color_image.make_general( XDF_F435W, ((850.0/435.0))*XDF_F850LP, ((1600.0/435.0))*XDF_F160W, os.path.join(outdir, 'XDF_RZH_Sept2013.png'), 3.0, 4.0,sigma_tuple=[0.01,0.01,0.01] )
    result = make_color_image.make_nasa( XDF_F435W_rotated, ((850.0/435.0))*XDF_F850LP_rotated, ((1600.0/435.0))*XDF_F160W_rotated, os.path.join(outdir, 'XDF_RZH_Dec2013_NASA.png'), 3.0, 4.0 )

    #5250*0.06 = 315 arcsec per side = 5.25 arcmin per side

    return 'made real HUDF image.'


def udfscript():






    #fwhm_pixels
    #

    #sigma

    hires_pixscale = 0.0416
    HST_pixscale = 0.06
    ratio = HST_pixscale/hires_pixscale
    sigma_tuple = ratio*np.asarray([0.0489, 0.0091, 0.0020])  #RMS noise per pixel estimated from Illingworth et al. 2013

    #turn off noise
    #sigma_tuple = sigma_tuple*0.0

    sigma_tuple = [0.01, 0.004, 0.0016]

    psfsig_arcsec = np.asarray([0.04,0.08,0.16])
    psfsig_pixels = psfsig_arcsec/hires_pixscale


    fwhm_pixels = 2.355*psfsig_pixels*0.0
    print 'psf fwhm, in pixels:  ', fwhm_pixels

    print sigma_tuple





    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER31_A"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC31_A.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER31_A1"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC31_A1.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER31_A2"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC31_A2.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER31_B"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC31_B.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER31_B1"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC31_B1.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER31_B2"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC31_B2.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    '''


    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER30_1e7_new_ism_both"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC30.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)





    print 'psf fwhm, in pixels:  ', fwhm_pixels

    print sigma_tuple

    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_zyx_COLDGAS/images/DECEMBER30_1e7_new_ism_both"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDC_DEC30.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)





    
    print 'psf fwhm, in pixels:  ', fwhm_pixels

    print sigma_tuple

    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_xyz_COLDGAS/images/DECEMBER30_1e7_new_ism_both"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDB_DEC30.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    



    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_zyx_COLDGAS/images/DECEMBER23_1e7_new_ism"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDC_DEC23_1e7_new_ism.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)





    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_zyx_COLDGAS/images/DECEMBER23_1e7_new_ism_both"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDC_DEC23_1e7_new_ism_both.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)





    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_zyx_COLDGAS/images/DECEMBER23_1e7_orig"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDC_DEC23_1e7_original.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)







    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_zyx_COLDGAS/images/DECEMBER23_3e7_new_ism_both"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDC_DEC23_3e7_new_ism_both.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)







    






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_xyz_COLDGAS/images/DECEMBER18_z00"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDB_DEC18_z00.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)
    



    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER18_z00"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC18_z00.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)







    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_zyx_COLDGAS/images/DECEMBER18_z00"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDC_DEC18_z00.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)








    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_xyz_COLDGAS/images/DECEMBER18_z05"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDB_DEC18_z05.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)
    



    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER18_z05"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC18_z05.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)



    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_zyx_COLDGAS/images/DECEMBER18_z05"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDC_DEC18_z05.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)










    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_xyz_COLDGAS/images/DECEMBER18_z03"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDB_DEC18_z03.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)
    



    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER18_z03"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC18_z03.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_zyx_COLDGAS/images/DECEMBER18_z03"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDC_DEC18_z03.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_xyz_COLDGAS/images/DECEMBER18_z01"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDB_DEC18_z01.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)
    



    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz_COLDGAS/images/DECEMBER18_z01"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDA_DEC18_z01.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_zyx_COLDGAS/images/DECEMBER18_z01"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[0].data  #1 = dandanxu,  2 = gfs, I think,  0 = no post-diffuse dust
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]    #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00_FIELDC_DEC18_z01.png',3.0,4.0,sigma_tuple=sigma_tuple,fwhm_pixels=fwhm_pixels)











    
    #SCI_BC03_MAPPINGS_ADAPTIVE_SEMIDUST

    ##2, 3, 6 ACS F435W, F606W, F850LP









    
    
    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n910FP/hudf_75Mpc_11_10/images/SCI_SB99_CF00_ADAPTIVE1"

    #read in data
    fitsfile = pyfits.open(os.path.join(folder,'bbdata.fits'))
    #z03_z12.info()
    imdata = fitsfile[0].data
    print imdata.shape

    filtertable = fitsfile[10]
    print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[3,:,:] ; r = imdata[6,:,:]  #F435W, F606W, F850LP  #filters_hstandjwst_new July 2013

    result = make_color_image.make_general(b,((606.0/435.0))*g,((850.0/435.0))*r,'SB99_CF00_ADAPTIVE_L75n910FP_ACSUDF.pdf',3.0,4.0,sigma_tuple=[0.01,0.01,0.01])
    '''


    
    return "Images made."



def shortscript():
    
    '''
    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz/images/"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'broadband_0056_cf00ish.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[12].data
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]  #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00ish_0056.png',3.0,4.0,sigma_tuple=[0.01,0.01,0.01])






    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz/images/"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'broadband_0056_cf00_test_withoutZtrend.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[12].data
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]  #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00ish_noZtrend_0056.png',3.0,4.0,sigma_tuple=[0.01,0.01,0.01])





    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz/images/"
    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'broadband_0056_cf00_test_withZtrend.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[12].data
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]  #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00ish_Ztrend_0056.png',3.0,4.0,sigma_tuple=[0.01,0.01,0.01])





    

    #read in data
    
    fitsfile = pyfits.open(os.path.join(folder,'broadband_0056_cf00_test_withoutZtrend_nodiffuse.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[12].data
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]  #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00ish_noZtrend_nodiffuse_0056.png',3.0,4.0,sigma_tuple=[0.01,0.01,0.01])





    
    fitsfile = pyfits.open(os.path.join(folder,'broadband_0056_cf00_test_1e8.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[12].data
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]  #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_general(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00ish_1e8_0056_LUPTON.png',3.0,4.0,sigma_tuple=[0.01,0.01,0.01])
    '''

    folder = "/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_yxz/images/"

    
    fitsfile = pyfits.open(os.path.join(folder,'broadband_nodiffuse_weakZtrend.fits'))

    #z03_z12.info()

    print fitsfile.info()
    imdata = fitsfile[12].data
    print imdata.shape
    print fitsfile[0].header.ascardlist()
    #print fitsfile[4].header.ascardlist()

    #filtertable = fitsfile[10]
    #print filtertable.data.field('filter')

    b = imdata[2,:,:] ; g = imdata[6,:,:] ; r = imdata[10,:,:]  #F435W, F850, F160W  #filters_hstandjwst_new July 2013

    result = make_color_image.make_nasa(b,((850.0/435.0))*g,((1600.0/435.0))*r,'BC03_CF00ish_1e8_nodiffuse_weakZtrend.png',3.0,4.0,sigma_tuple=[0.01,0.01,0.01])



    return "Short Image Made."




if __name__=="__main__":

    result = udfscript()
    
    print result



    #result = shortscript()

    #print result

    #result = realdatascript()

    #print result
