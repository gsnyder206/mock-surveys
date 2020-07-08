

import math
import astropy
import astropy.io.fits as fits
import matplotlib
import matplotlib.pyplot as pyplot
import os
import scipy
import scipy.ndimage
import scipy.ndimage.interpolation
import numpy as np

import congrid #http://scipy-cookbook.readthedocs.io/items/Rebinning.html#Example-3
import make_color_image #written by Gregory Snyder; Also in gsnyder206/synthetic-image-morph

sq_arcsec_per_sr = 42545170296.0
c = 3.0e8

hlsp_dir='/astro/snyder_lab2/New_HydroART_images/VELA_v2/luvoir_mocks/luvoir_mosaics/'
#/astro/snyder_lab2/New_HydroART_images/VELA_v2/luvoir_mocks/luvoir_mosaics/hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f435w_v1_lightcone.fits


def do_hst(fieldstr,alph,Q,rf,gf,bf):
    pass


def load_luvoir(fieldstr,subim='00',smooth=''):
    
    #in nJy
    fielda_f336 =fits.open(hlsp_dir+'hlsp_misty_mockluvoir_imager_'+fieldstr+'-subimage'+subim+'-9-8'+smooth+'_f336w_v1_lightcone.fits')[0].data
    fielda_f435 =fits.open(hlsp_dir+'hlsp_misty_mockluvoir_imager_'+fieldstr+'-subimage'+subim+'-9-8'+smooth+'_f435w_v1_lightcone.fits')[0].data
    fielda_f606 =fits.open(hlsp_dir+'hlsp_misty_mockluvoir_imager_'+fieldstr+'-subimage'+subim+'-9-8'+smooth+'_f606w_v1_lightcone.fits')[0].data

    #f775w corrupted for some reason?  many nans
    #fielda_f775 =fits.open(hlsp_dir+'hlsp_misty_mockluvoir_imager_'+fieldstr+'-subimage'+subim+'-9-8'+smooth+'_f775w_v1_lightcone.fits')[0].data
    fielda_f850 =fits.open(hlsp_dir+'hlsp_misty_mockluvoir_imager_'+fieldstr+'-subimage'+subim+'-9-8'+smooth+'_f850lp_v1_lightcone.fits')[0].data

    fielda_f115 =fits.open(hlsp_dir+'hlsp_misty_mockluvoir_imager_'+fieldstr+'-subimage'+subim+'-9-8'+smooth+'_f115w_v1_lightcone.fits')[0].data
    fielda_f150 =fits.open(hlsp_dir+'hlsp_misty_mockluvoir_imager_'+fieldstr+'-subimage'+subim+'-9-8'+smooth+'_f150w_v1_lightcone.fits')[0].data
    fielda_f200 =fits.open(hlsp_dir+'hlsp_misty_mockluvoir_imager_'+fieldstr+'-subimage'+subim+'-9-8'+smooth+'_f200w_v1_lightcone.fits')[0].data

    pixsize=fielda_f200 =fits.open(hlsp_dir+'hlsp_misty_mockluvoir_imager_'+fieldstr+'-subimage'+subim+'-9-8'+smooth+'_f200w_v1_lightcone.fits')[0].header['PIXSIZE']

    ra= rf * (0.333*fielda_f115 + 0.333*fielda_f150 + 0.333*fielda_f200)
    ga= gf * (0.50*fielda_f850 + 0.5*fielda_f606)
    ba= bf * (0.50*fielda_f435 + 0.50*fielda_f336)

    return ra,ga,ba,pixsize

def do_luvoir(fieldstr,alph,Q,rf,gf,bf,subim='00',smooth='',diameter=None):

    print('loading data... ')
    
    ra,ga,ba,pix_arcsec=load_luvoir(fieldstr,subim=subim,smooth=smooth)

    print('loaded data, making rgb... ')

    
    if diameter is None:
        rgb_field = make_color_image.make_interactive_nasa(ba,ga,ra,alph,Q)
        psfstr=''
    else:
        fwhm_radians=1.22*np.asarray([0.40,0.75,1.5])*(1.0e-6)/diameter
        fwhm_arcsec=fwhm_radians*(180.0/math.pi)*3600.0
        print('applying PSF FWHM, in arcsec: ', fwhm_arcsec)
        fwhm_pixels=fwhm_arcsec/pix_arcsec
        rgb_field= make_color_image.make_interactive_nasa(ba,ga,ra,alph,Q,fwhm_pixels=fwhm_pixels)
        psfstr='_'+str(diameter)+'meter'
        
    print('created color image, saving...  ')
    
    f1 = pyplot.figure(figsize=(8.192,8.192), dpi=1000)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)

    axi=f1.add_subplot(111)
    axi.imshow(rgb_field,interpolation='nearest',aspect='auto',origin='lower')

    f1.savefig('luvoir_render_'+fieldstr+'-subimage'+subim+smooth+psfstr+'.png',dpi=1000)
    pyplot.close(f1)



if __name__=="__main__":
    rf= 1.3
    gf= 0.8
    bf= 1.3
    alph=10.0
    Q=2.0
    


    do_luvoir('fielda',alph,Q,rf,gf,bf,subim='00')
    do_luvoir('fielda',alph,Q,rf,gf,bf,subim='00',smooth='-smooth')
    
    do_luvoir('fielda',alph,Q,rf,gf,bf,subim='00',smooth='-smooth',diameter=15.0)
    do_luvoir('fielda',alph,Q,rf,gf,bf,subim='00',smooth='-smooth',diameter=2.5)

    do_luvoir('fielda',alph,Q,rf,gf,bf,subim='00',diameter=15.0)
    do_luvoir('fielda',alph,Q,rf,gf,bf,subim='00',diameter=2.5)


    #do_luvoir('fielda',alph,Q,rf,gf,bf,subim='00',smooth='-smooth',diameter=6.5)

    #do_jwst('fielda',alph,Q,rf,gf,bf)
    
    

