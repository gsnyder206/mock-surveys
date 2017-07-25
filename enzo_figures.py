import illustris_python as ilpy
import h5py
import math
import string
import sys
import struct
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.colors as pycolors
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import numpy.ma as ma
import asciitable
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import scipy.odr as odr
from scipy.integrate import quad
import glob
import os
import gzip
import shutil
import astropy.io.ascii as ascii
import warnings
import subprocess
import astropy
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.cosmology import WMAP9,z_at_value
import copy
import datetime
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import PolyCollection
import pandas
import make_color_image
import random
import gfs_twod_histograms as gth

class snapshot:
    def __init__(self):
        return



def jwst_nirspec_figure(alph=0.1,Q=1.0,dx=0,dy=0,Npix=400.0,do_psf=False):

    natch_z2='nref11/RD0020/mcrx.fits'
    extra_z2='nref11_refine200kpc_z4to2/RD0020/mcrx.fits'
    extra_z3='nref11_refine200kpc_z4to2/RD0017/mcrx.fits'

    cs='CAMERA3-NONSCATTER'

    natch_z2_o=snapshot()
    extra_z2_o=snapshot()
    extra_z3_o=snapshot()
    

    rf=21 #F200W         0.75
    gf=19 #F115W         0.40
    bf=4  #ACS/F775W   rest far-UV   0.25

    il=1310
    in2=1316
    ic=1296
    
    fudge=1.0
    lfact=1.0

    mid = np.int64(400.0/2)
    if Npix is not None:
        delt=np.int64(Npix/2)
    

    filen='nirspec.pdf'

    
    extra_z2_o.ha =fudge*pyfits.open(extra_z2)[cs].data[il,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*lfact
    extra_z3_o.ha =fudge*pyfits.open(extra_z3)[cs].data[il,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*lfact

    extra_z2_o.n2 =fudge*pyfits.open(extra_z2)[cs].data[in2,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*lfact
    extra_z3_o.n2 =fudge*pyfits.open(extra_z3)[cs].data[in2,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*lfact

    extra_z2_o.c =fudge*pyfits.open(extra_z2)[cs].data[ic,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*lfact
    extra_z3_o.c =fudge*pyfits.open(extra_z3)[cs].data[ic,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*lfact

    
    natch_z2_o.ha =fudge*pyfits.open(natch_z2)[cs].data[il,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*lfact

    natch_z2_o.n2 =fudge*pyfits.open(natch_z2)[cs].data[in2,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*lfact

    natch_z2_o.c =fudge*pyfits.open(natch_z2)[cs].data[ic,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*lfact

    

    kpc_per_arcsec_z2=WMAP9.kpc_proper_per_arcmin(2.0).value/60.0
    kpc_per_arcsec_z3=WMAP9.kpc_proper_per_arcmin(2.75).value/60.0
    kpc_per_pix=0.125
    
    

    #high-res version first

    
    spinelw=1.0
    
    f1 = pyplot.figure(figsize=(10.0,10.0), dpi=200)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)

    axi=f1.add_subplot(2,2,1)
    axi.set_xticks([]) ; axi.set_yticks([])

    axi.annotate('JWST-Nirspec/IFU',(0.50,0.88),xycoords='axes fraction',color='Black',ha='center',va='center',fontsize=25)

    rat=np.log10(natch_z2_o.n2/natch_z2_o.ha)


    vmin=-0.5
    vmax=0.0
    themap = cm.viridis

    rat=np.where(rat > vmax, -1000.0, rat)
    Zm = ma.masked_where( rat < -100.0, rat)
    cscale_function = pycolors.Normalize(vmin=vmin,vmax=vmax,clip=True)
    
    
    axi.imshow(np.fliplr(np.transpose(Zm)),interpolation='nearest',origin='upper',vmin=vmin,vmax=vmax,cmap=themap)
    

    for ss in axi.spines:
        s=axi.spines[ss]
        s.set_color('white')
        s.set_linewidth(spinelw)


        
    axi=f1.add_subplot(2,2,2)

    axi.set_xticks([]) ; axi.set_yticks([])

    rat=np.log10(extra_z2_o.n2/extra_z2_o.ha)


    vmin=-0.5
    vmax=0.0
    themap = cm.viridis

    rat=np.where(rat > vmax, -1000.0, rat)
    Zm = ma.masked_where( rat < -100.0, rat)
    cscale_function = pycolors.Normalize(vmin=vmin,vmax=vmax,clip=True)
    
    
    axi.imshow(np.fliplr(np.transpose(Zm)),interpolation='nearest',origin='upper',vmin=vmin,vmax=vmax,cmap=themap)

    
    for ss in axi.spines:
        s=axi.spines[ss]
        s.set_color('white')
        s.set_linewidth(spinelw)




        
    axi=f1.add_subplot(2,2,4)
    axi.set_xticks([]) ; axi.set_yticks([])

    rat=np.log10(extra_z3_o.n2/extra_z3_o.ha)

    vmin=-0.5
    vmax=0.0
    themap = cm.viridis

    rat=np.where(rat > vmax, -1000.0, rat)
    Zm = ma.masked_where( rat < -100.0, rat)
    cscale_function = pycolors.Normalize(vmin=vmin,vmax=vmax,clip=True)
    
    
    obj=axi.imshow(np.fliplr(np.transpose(Zm)),interpolation='nearest',origin='upper',vmin=vmin,vmax=vmax,cmap=themap)

    gth.make_colorbar(obj,f1,title=r'$log_{10} [NII]/H \alpha $',ticks=[-0.5,0.0],loc=[0.60,0.50,0.3,0.05],fontsize=25)

    
    for ss in axi.spines:
        s=axi.spines[ss]
        s.set_color('white')
        s.set_linewidth(spinelw)






    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)


        
    f1.savefig(filen,dpi=200)
    pyplot.close(f1)
    

    return




    

def jwst_image_figure(alph=0.1,Q=1.0,dx=0,dy=0,Npix=1000.0,do_psf=False,cs='CAMERA5-BROADBAND-NONSCATTER',label=''):


    nref11_z2='nref11/RD0020/broadbandz.fits'
    nref11_z3='nref11/RD0017/broadbandz.fits'

    extra_z2='nref11_refine200kpc_z4to2/RD0020/broadbandz.fits'
    extra_z3='nref11_refine200kpc_z4to2/RD0017/broadbandz.fits'


    nref11_z2_o=snapshot()
    nref11_z3_o=snapshot()
    extra_z2_o=snapshot()
    extra_z3_o=snapshot()
    

    rf=21 #F200W         0.75
    gf=19 #F115W         0.40
    bf=4  #ACS/F775W   rest far-UV   0.25

    fudge=1.0

    bfact=1.5
    gfact=0.5/(0.115**2)
    rfact=2.0/(0.200**2)

    mid = np.int64(1000.0/2)
    if Npix is not None:
        delt=np.int64(Npix/2)
    
    nref11_z2_o.r=fudge*pyfits.open(nref11_z2)[cs].data[rf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*rfact
    nref11_z2_o.g=fudge*pyfits.open(nref11_z2)[cs].data[gf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*gfact
    nref11_z2_o.b=fudge*pyfits.open(nref11_z2)[cs].data[bf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*bfact

    nref11_z3_o.r=fudge*pyfits.open(nref11_z3)[cs].data[rf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*rfact
    nref11_z3_o.g=fudge*pyfits.open(nref11_z3)[cs].data[gf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*gfact
    nref11_z3_o.b=fudge*pyfits.open(nref11_z3)[cs].data[bf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*bfact


    print(np.max(nref11_z2_o.r),np.max(nref11_z2_o.b))
    
    extra_z2_o.r=fudge*pyfits.open(extra_z2)[cs].data[rf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*rfact
    extra_z2_o.g=fudge*pyfits.open(extra_z2)[cs].data[gf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*gfact
    extra_z2_o.b=fudge*pyfits.open(extra_z2)[cs].data[bf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*bfact

    extra_z3_o.r=fudge*pyfits.open(extra_z3)[cs].data[rf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*rfact
    extra_z3_o.g=fudge*pyfits.open(extra_z3)[cs].data[gf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*gfact
    extra_z3_o.b=fudge*pyfits.open(extra_z3)[cs].data[bf,dx+mid-delt:dx+mid+delt,dy+mid-delt:dy+mid+delt]*bfact


    
    print(extra_z2_o.r.shape)

    kpc_per_arcsec_z2=WMAP9.kpc_proper_per_arcmin(2.0).value/60.0
    kpc_per_arcsec_z3=WMAP9.kpc_proper_per_arcmin(2.75).value/60.0
    kpc_per_pix=0.05
    
    
    if do_psf==False:
        filen='images_hires_'+label+'_'+cs[0:7]+'.pdf'
        nref11_z2_o.rgbthing = make_color_image.make_interactive(nref11_z2_o.b,nref11_z2_o.g,nref11_z2_o.r,alph,Q)
        nref11_z3_o.rgbthing = make_color_image.make_interactive(nref11_z3_o.b,nref11_z3_o.g,nref11_z3_o.r,alph,Q)
        extra_z2_o.rgbthing = make_color_image.make_interactive(extra_z2_o.b,extra_z2_o.g,extra_z2_o.r,alph,Q)
        extra_z3_o.rgbthing = make_color_image.make_interactive(extra_z3_o.b,extra_z3_o.g,extra_z3_o.r,alph,Q)
    else:
        filen='images_psf_'+label+'_'+cs[0:7]+'.pdf'
        psf_fwhm_arcsec=0.05
        fwhm_pixels_z2=psf_fwhm_arcsec*kpc_per_arcsec_z2/kpc_per_pix
        fwhm_pixels_z3=psf_fwhm_arcsec*kpc_per_arcsec_z3/kpc_per_pix
        print('KPC per as, z2: ', kpc_per_arcsec_z2)
        print('FWHM arcsec:    ', psf_fwhm_arcsec)
        print('FWHM pixels z2: ', fwhm_pixels_z2)
        print('FWHM pixels z3: ', fwhm_pixels_z3)
       
        nref11_z2_o.rgbthing = make_color_image.make_interactive(nref11_z2_o.b,nref11_z2_o.g,nref11_z2_o.r,alph,Q,fwhm_pixels=fwhm_pixels_z2)
        nref11_z3_o.rgbthing = make_color_image.make_interactive(nref11_z3_o.b,nref11_z3_o.g,nref11_z3_o.r,alph,Q,fwhm_pixels=fwhm_pixels_z3)
        extra_z2_o.rgbthing = make_color_image.make_interactive(extra_z2_o.b,extra_z2_o.g,extra_z2_o.r,alph,Q,fwhm_pixels=fwhm_pixels_z2)
        extra_z3_o.rgbthing = make_color_image.make_interactive(extra_z3_o.b,extra_z3_o.g,extra_z3_o.r,alph,Q,fwhm_pixels=fwhm_pixels_z3)        
    
    #high-res version first

    
    spinelw=1.0
    
    f1 = pyplot.figure(figsize=(10.0,10.0), dpi=200)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)

    
    axi=f1.add_subplot(2,2,1)
    axi.set_xticks([]) ; axi.set_yticks([])
    axi.imshow(nref11_z2_o.rgbthing,interpolation='nearest',origin='lower')

    axi.annotate('HST-ACS/F775W',(0.25,0.94),xycoords='axes fraction',color='Blue',ha='center',va='center',fontsize=20)
    axi.annotate('JWST-NC/F115W',(0.25,0.88),xycoords='axes fraction',color='Green',ha='center',va='center',fontsize=20)
    axi.annotate('JWST-NC/F200W',(0.25,0.82),xycoords='axes fraction',color='Red',ha='center',va='center',fontsize=20)

    #axi.annotate(anfk,(0.5,0.15),xycoords='axes fraction',color=anc,ha='center',va='center',fontsize=12,backgroundcolor='Black')

    for ss in axi.spines:
        s=axi.spines[ss]
        s.set_color('white')
        s.set_linewidth(spinelw)


        
    axi=f1.add_subplot(2,2,2)

    axi.set_xticks([]) ; axi.set_yticks([])
    axi.imshow(extra_z2_o.rgbthing,interpolation='nearest',origin='lower')

    #10kpc scale bar
    if label=='zoom':
        axi.plot([delt-10,delt+10],[delt/3,delt/3],marker='None',linestyle='solid',color='White',lw=5)
        axi.annotate('1 kpc', (delt,delt/4), xycoords='data',color='White',ha='center',va='center',fontsize=20)
        axi.annotate('z = 2', (0.5*delt,2*delt*0.90), xycoords='data',color='White',ha='center',va='center',fontsize=20)
    else:
        axi.plot([delt-100,delt+100],[delt/3,delt/3],marker='None',linestyle='solid',color='White',lw=5)
        axi.annotate('10 kpc', (delt,delt/4), xycoords='data',color='White',ha='center',va='center',fontsize=20)
        axi.annotate('z = 2', (0.5*delt,2*delt*0.90), xycoords='data',color='White',ha='center',va='center',fontsize=20)
    
    if do_psf==False:
        axi.annotate('No PSF', (1.5*delt,2*delt*0.90), xycoords='data',color='White',ha='center',va='center',fontsize=20)
    else:
        axi.annotate('0.05" FWHM', (1.5*delt,2*delt*0.90), xycoords='data',color='White',ha='center',va='center',fontsize=20)
        
    for ss in axi.spines:
        s=axi.spines[ss]
        s.set_color('white')
        s.set_linewidth(spinelw)


        
    axi=f1.add_subplot(2,2,3)
    axi.imshow(nref11_z3_o.rgbthing,interpolation='nearest',origin='lower')

    axi.set_xticks([]) ; axi.set_yticks([])
    axi.annotate('natural refine', (delt,0.25*delt), xycoords='data',color='White',ha='center',va='center',fontsize=20)

    for ss in axi.spines:
        s=axi.spines[ss]
        s.set_color('white')
        s.set_linewidth(spinelw)




        
    axi=f1.add_subplot(2,2,4)
    axi.set_xticks([]) ; axi.set_yticks([])
    axi.imshow(extra_z3_o.rgbthing,interpolation='nearest',origin='lower')

    axi.annotate('z = 2.75', (0.5*delt,2*delt*0.90), xycoords='data',color='White',ha='center',va='center',fontsize=20)
    axi.annotate('forced refine', (delt,0.25*delt), xycoords='data',color='White',ha='center',va='center',fontsize=20)

    
    for ss in axi.spines:
        s=axi.spines[ss]
        s.set_color('white')
        s.set_linewidth(spinelw)






    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)


        
    f1.savefig(filen,dpi=200)
    pyplot.close(f1)
    

    return






if __name__=="__main__":
    jwst_image_figure()
    jwst_nirspec_figure()
