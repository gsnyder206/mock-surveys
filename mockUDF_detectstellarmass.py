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

sq_arcsec_per_sr = 42545170296.0
c = 3.0e8



if __name__=="__main__":


    im_directory = '/Users/gsnyder/Documents/Projects/mockUDFs/FITS_FILES/fromOdyssey/YXZ/DECEMBER30_1e7_new_ism_both'
    bbdata = pyfits.open(os.path.join(im_directory,'bbdata.fits'))

    mstar_grid = (bbdata['STELLAR_MASS'].data)
    print mstar_grid.shape

    se_dir = '/Users/gsnyder/Documents/Projects/mockUDFs/SE_RUNS/FIELDA'
    se_runs = ['FIELDA_SB25','FIELDA_SB27','FIELDA_SB29']
    mag_lims = [25.0,27.0,29.0]


    imfilename = 'WideFieldDepthGrid.pdf'
    f9 = pyplot.figure(figsize=(12.0,8.0), dpi=1000)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)

    FWHM_arcsec = 0.115
    Hl = 1.60
    Jl = 1.25
    Zl = 0.900
    Il = 0.775
    Vl = 0.606
    lambdas = np.asarray([Hl,Jl,Zl,Il,Vl])

    r_i = 0 ; g_i = 1 ; b_i = 2
    gscale = 0.85*lambdas[b_i]/lambdas[g_i]
    rscale = lambdas[b_i]/lambdas[r_i]
    deltax = 1.0/3.0
    deltay = 1.0/2.0

    factor = 2.814
    cuti = 0#int(4096.0/3.0)
    cutj = cuti
    dij = int(4096.0/factor)

    i=0
    for se in se_runs:
        #print se
        startx = 0.0 + i*0.33333333
        starty = 0.50
        starty2 = 0.0

        se_folder = os.path.join(se_dir,se)
        segmapfile = os.path.join(se_folder,'FIELDA_hphotom_comb_seg.fits')
        hmapfile = os.path.join(se_folder,se+'_F160W.fits')
        jmapfile = os.path.join(se_folder,se+'_F125W.fits')
        zmapfile = os.path.join(se_folder,se+'_F850LP.fits')

        segmap = pyfits.open(segmapfile)[0].data
        nzi = np.nonzero(segmap)
        #nzi2 = np.where(segmap > 1.0e-5)[0]
        #segmap2 = segmap*1.0
        #segmap2[nzi2] = np.ones_like(segmap2[nzi2])

        #segmap[nzi] = np.ones_like(segmap[nzi])

        hmap = pyfits.open(hmapfile)[0].data
        jmap = pyfits.open(jmapfile)[0].data
        zmap = pyfits.open(zmapfile)[0].data
        header = pyfits.open(hmapfile)[0].header

        target_ratio = 10.0**(-0.4*(28.0-mag_lims[i]))
        fluxscale = target_ratio*3.0e-14
        b = zmap*fluxscale
        g = jmap*fluxscale
        r = hmap*fluxscale

        alph = 5.0
        Q = 9.0

        sigma_nJySr = header['RMS']
        pixel_Sr = (header['PIXSCALE']**2)/sq_arcsec_per_sr
        sigma_nJy = sigma_nJySr*pixel_Sr
        N = (3.0*FWHM_arcsec/header['PIXSCALE'])**2
        mag_equiv = -2.5*np.log10(5.0*sigma_nJy*(1.0e-9)/(3631.0/(N**0.5)))
        mag_string = '{:4.1f}'.format(mag_equiv)

        slw = 2.0

        rgbdata = make_color_image.make_interactive_nasa(b,g*gscale,r*rscale,alph,Q)
        axi = pyplot.axes([startx,starty,deltax,deltay],frameon=True,axisbg='black')
        axi.set_xticks([]) ; axi.set_yticks([])
        axi.imshow(rgbdata,interpolation='nearest',origin='lower',extent=[0,1,0,1])  
        axi.annotate(mag_string,(0.5,0.09),xycoords='axes fraction',ha='center',va='center',color='white',size=12)

        axi.add_patch(pyplot.Rectangle((0.01,0.01),1.0/factor,1.0/factor,facecolor='None',edgecolor='LightGreen',linewidth=slw,linestyle='dotted'))
        axi.annotate('1 arcmin',(0.5/factor,1.0/factor+0.05),xycoords='axes fraction',ha='center',va='center',color='LightGreen',size=12)

        rgbdata2 = make_color_image.make_interactive_nasa(b[cuti:cuti+dij,cuti:cuti+dij],g[cuti:cuti+dij,cuti:cuti+dij]*gscale,r[cuti:cuti+dij,cuti:cuti+dij]*rscale,alph,Q)
        axi2 = pyplot.axes([startx,starty2,deltax,deltay],frameon=True,axisbg='black')
        axi2.set_xticks([]) ; axi2.set_yticks([])
        axi2.imshow(rgbdata2,interpolation='nearest',origin='lower',extent=[-1,1,-1,1])  
        for pos in ['top', 'bottom', 'right', 'left']:
            axi2.spines[pos].set_edgecolor('LightGreen')
            axi2.spines[pos].set_linewidth(slw)


        mstar_detected = np.sum( mstar_grid[nzi])
        hflux_detected = np.sum( hmap[nzi])
        zflux_detected = np.sum( zmap[nzi])

        totalh = np.sum(hmap)
        totalz = np.sum(zmap)
        #totalsm = np.sum(mstar_grid)

        #segmap2 = segmap*1.0
        #segmap2[nzi] = np.ones_like(segmap[nzi])


        #print segmap2.shape
        #newHDU = pyfits.PrimaryHDU(segmap2)
        #newlist = pyfits.HDUList([newHDU])
        #newlist.writeto(se+'_allsegmap.fits',clobber=True)

        print se, mstar_detected, hflux_detected, totalh, zflux_detected, totalz #np.sum(mstar_grid*segmap2), np.sum(hmap*segmap2)
        i=i+1
        pyplot.annotate('2 < z < 3',(0.05,0.97),xycoords='figure fraction',ha='center',va='center',color='white',size=16,backgroundcolor='black')


    f9.savefig(imfilename,format='pdf',dpi=1000)
    pyplot.close(f9)
