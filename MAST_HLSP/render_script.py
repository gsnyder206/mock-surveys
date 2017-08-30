


import astropy
import astropy.io.fits as fits
import matplotlib
import matplotlib.pyplot as pyplot
import os
import scipy
import scipy.ndimage
import scipy.ndimage.interpolation

import congrid #http://scipy-cookbook.readthedocs.io/items/Rebinning.html#Example-3
import make_color_image #written by Gregory Snyder; Also in gsnyder206/synthetic-image-morph

sq_arcsec_per_sr = 42545170296.0
c = 3.0e8

hlsp_dir='/astro/snyder_lab2/Illustris/Lightcones/Lightcone_Catalog_Images/'



def do_hst_illustris(fieldstr,alph,Q,rf,gf,bf):
    
    #in nJy
    fielda_f435 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-acs_f435w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')['IMAGE_PSF'].data
    fielda_f606 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-acs_f606w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')['IMAGE_PSF'].data
    fielda_f775 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-acs_f775w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')['IMAGE_PSF'].data
    fielda_f814 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-acs_f814w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')['IMAGE_PSF'].data
    fielda_f850 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-acs_f850lp_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')['IMAGE_PSF'].data
    
    fielda_f105 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-wfc3_f105w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')['IMAGE_PSF'].data
    fielda_f125 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-wfc3_f125w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')['IMAGE_PSF'].data
    fielda_f140 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-wfc3_f140w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')['IMAGE_PSF'].data
    fielda_f160 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-wfc3_f160w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')['IMAGE_PSF'].data

    ra= rf * (0.25*fielda_f105 +0.25*fielda_f125 + 0.25*fielda_f140 + 0.25*fielda_f160)
    ga= gf * (0.50*fielda_f850 + 0.50*fielda_f775)
    ba= bf * (0.50*fielda_f435 + 0.50*fielda_f606)

    gnew=congrid.congrid(ga,ra.shape)*4.0  #preserve surface brightness
    bnew=congrid.congrid(ba,ra.shape)*4.0

    print(fielda_f435.shape)
    print(ra.shape,ga.shape,gnew.shape)
    
    rgb_field = make_color_image.make_interactive_nasa(bnew,gnew,ra,alph,Q)


    f1 = pyplot.figure(figsize=(10.0,10.0), dpi=600)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)

    axi=f1.add_subplot(111)
    axi.imshow(rgb_field,interpolation='nearest',aspect='auto',origin='lower')

    f1.savefig('illustris_render_'+fieldstr+'.pdf',dpi=600)
    pyplot.close(f1)


def do_jwst_illustris(fieldstr,alph,Q,rf,gf,bf,x=None,y=None,n=None):
    
    #in nJy
    ui='IMAGE_PSF'
    
    fielda_f435 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-acs_f435w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data
    fielda_f606 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-acs_f606w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data
    fielda_f775 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_hst-acs_f775w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data
    
    fielda_f090 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_jwst-nircam_f090w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data

    fielda_f115 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_jwst-nircam_f115w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data
    fielda_f150 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_jwst-nircam_f150w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data
    fielda_f200 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_jwst-nircam_f200w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data
    
    fielda_f277 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_jwst-nircam_f277w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data
    fielda_f356 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_jwst-nircam_f356w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data
    fielda_f444 =fits.open(hlsp_dir+'mag30-'+fieldstr+'-11-10/hlsp_misty_illustris_jwst-nircam_f444w_mag30-'+fieldstr+'-11-10_v1_lightcone.fits')[ui].data

    ra= rf * (0.50*fielda_f150 + 0.50*fielda_f200)
    ga= gf * (0.50*fielda_f090 + 0.50*fielda_f115)
    ba= bf * (0.33*fielda_f435 + 0.33*fielda_f606 + 0.33*fielda_f775)

    gnew=congrid.congrid(ga,ra.shape)*1.0  #preserve surface brightness
    bnew=congrid.congrid(ba,ra.shape)*1.0

    print(fielda_f435.shape)
    print(ra.shape,ga.shape,gnew.shape)


    if n is not None:
        rgb_field = make_color_image.make_interactive_nasa(bnew[x:x+n,y:y+n],gnew[x:x+n,y:y+n],ra[x:x+n,y:y+n],alph,Q)
    else:
        rgb_field = make_color_image.make_interactive_nasa(bnew,gnew,ra,alph,Q)


    f1 = pyplot.figure(figsize=(10.0,10.0), dpi=600)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)

    axi=f1.add_subplot(111)
    axi.imshow(rgb_field,interpolation='nearest',aspect='auto',origin='lower')

    if n is not None:
        saven='illustris_jwstzoom_'+fieldstr+'.pdf'
    else:
        saven='illustris_jwst_'+fieldstr+'.pdf'
    
    f1.savefig(saven,dpi=600)
    pyplot.close(f1)
    



def do_xdf(alph,Q,rf,gf,bf):

    xdf_file='/astro/snyder_lab2/Illustris/Illustris_MockUDFs/L75n1820FP/XDF/HST_XDF2012_SunriseUnits_60mas.fits'

    xdf_hdus=fits.open(xdf_file)
    print(xdf_hdus.info())

    pixsize_arcsec=0.06
    pixel_Sr = (pixsize_arcsec**2)/sq_arcsec_per_sr 
    to_nJy_per_Sr_perlambdasqum = (1.0e9)*(1.0e14)/c   
    to_nJy_per_pix = to_nJy_per_Sr_perlambdasqum*pixel_Sr

    dx=200 ; dy=200
    xx1=2200+dx 
    xx2=5042+dx
    xy1=2200+dy
    xy2=5042+dy
    #gives 2.841 arcmin
    
    xdf_f435 = scipy.ndimage.interpolation.rotate( fits.open(xdf_file)['F435W'].data, 42)[xx1:xx2,xy1:xy2] * to_nJy_per_pix * (0.435)**2
    xdf_f606 = scipy.ndimage.interpolation.rotate( fits.open(xdf_file)['F606W'].data, 42)[xx1:xx2,xy1:xy2] * to_nJy_per_pix * (0.606)**2
    xdf_f775 = scipy.ndimage.interpolation.rotate( fits.open(xdf_file)['F775W'].data, 42)[xx1:xx2,xy1:xy2] * to_nJy_per_pix * (0.775)**2
    xdf_f814 = scipy.ndimage.interpolation.rotate( fits.open(xdf_file)['F814W'].data, 42)[xx1:xx2,xy1:xy2] * to_nJy_per_pix * (0.814)**2
    xdf_f850 = scipy.ndimage.interpolation.rotate( fits.open(xdf_file)['F850LP'].data, 42)[xx1:xx2,xy1:xy2] * to_nJy_per_pix * (0.850)**2

    xdf_f105 = scipy.ndimage.interpolation.rotate( fits.open(xdf_file)['F105W'].data, 42)[xx1:xx2,xy1:xy2] * to_nJy_per_pix * (1.05)**2
    xdf_f125 = scipy.ndimage.interpolation.rotate( fits.open(xdf_file)['F125W'].data, 42)[xx1:xx2,xy1:xy2] * to_nJy_per_pix * (1.25)**2
    xdf_f140 = scipy.ndimage.interpolation.rotate( fits.open(xdf_file)['F140W'].data, 42)[xx1:xx2,xy1:xy2] * to_nJy_per_pix * (1.40)**2
    xdf_f160 = scipy.ndimage.interpolation.rotate( fits.open(xdf_file)['F160W'].data, 42)[xx1:xx2,xy1:xy2] * to_nJy_per_pix * (1.60)**2

    print(xdf_f435.shape)
    

    rx= rf * (0.25*xdf_f105 +0.25*xdf_f125 + 0.25*xdf_f140 + 0.25*xdf_f160)
    gx= gf * (0.50*xdf_f850 + 0.50*xdf_f775)
    bx= bf * (0.50*xdf_f435 + 0.50*xdf_f606)



    rgb_xdf = make_color_image.make_interactive_nasa(bx,gx,rx,alph,Q)


    f1 = pyplot.figure(figsize=(10.0,10.0), dpi=600)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)

    axi=f1.add_subplot(111)
    axi.imshow(rgb_xdf,interpolation='nearest',aspect='auto',origin='lower')

    f1.savefig('xdf_render.pdf',dpi=600)
    pyplot.close(f1)
    

    return



if __name__=="__main__":
    rf= 1.0
    gf= 1.0
    bf= 1.0
    alph=1.0
    Q=6.0
    
    #plot real XDF
    #do_xdf(alph,Q,rf,gf,bf)



    #Illustris XDF analogue

    #do_hst_illustris('fielda',alph,Q,rf,gf,bf)
    #do_hst_illustris('fieldb',alph,Q,rf,gf,bf)
    #do_hst_illustris('fieldc',alph,Q,rf,gf,bf)



    
    
    #Illustris HST+ JWST

    alph=2.0
    Q=7.0
    rf=1.3 ; gf=1 ; bf=1.2
    
    #do_jwst_illustris('fielda',alph,Q,rf,gf,bf)
    #do_jwst_illustris('fieldb',alph,Q,rf,gf,bf)
    #do_jwst_illustris('fieldc',alph,Q,rf,gf,bf)

    do_jwst_illustris('fieldc',alph,Q,rf,gf,bf,x=1200,y=4600,n=600)
    
    

