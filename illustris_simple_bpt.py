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
from astropy.cosmology import WMAP7,z_at_value
import copy
import datetime
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import PolyCollection
import PyML
import PyML.machinelearning
import gfs_twod_histograms as gth
import pandas
import make_color_image
import photutils
import gfs_sublink_utils as gsu
import random

ilh = 0.704
illcos = astropy.cosmology.FlatLambdaCDM(H0=70.4,Om0=0.2726,Ob0=0.0456)

shfields=['SubhaloBHMass','SubhaloBHMdot','SubhaloMass','SubhaloMassInRadType','SubhaloMassType','SubhaloSFR','SubhaloSFRinRad','SubhaloStarMetallicity','SubhaloGasMetallicitySfrWeighted']

def get_line_fluxes(sf='cont',n='n10',Z='1.0',q='8.0e7',M='K'):
    filen = 'SB99_'+sf+'_'+n+'/spec_Z'+Z+'_8Myr_q'+q+'_SB99_'+M+'.ph4.gz'
    dat = ascii.read(filen,data_start=20)

    io3=np.where(dat['col1']==5006.77)[0][0]
    iha=np.where(dat['col1']==6562.8)[0][0]
    in2=np.where(dat['col1']==6583.34)[0][0]


    o3hb = dat['col3'][io3] ; print(dat['col1'][io3])
    hahb = dat['col3'][iha] ; print(dat['col1'][iha])
    n2hb = dat['col3'][in2] ; print(dat['col1'][in2])

    return o3hb,hahb,n2hb

def get_subhalos(snapnum=75,mstarlim=10.0**10.5):
    #mstar_1 = dataobjs['subhalos'][s]['SubhaloMassInRadType'][:,4].flatten()*(1.0e10)/ilh
    #sfr_1 = dataobjs['subhalos'][s]['SubhaloSFR'][:].flatten()
    #mhalo_1 = dataobjs['subhalos'][s]['SubhaloMass'][:].flatten()*(1.0e10)/ilh
    #bhrate_1 = dataobjs['subhalos'][s]['SubhaloBHMdot'][:].flatten()*((1.0e10)/ilh)/(0.978*1.0e9/ilh)
    #bhmass_1 = dataobjs['subhalos'][s]['SubhaloBHMass'][:].flatten()*(1.0e10)/ilh 

    subhalo_dict = ilpy.groupcat.loadSubhalos('/astro/snyder_lab2/Illustris/Illustris-1',snapnum,fields=shfields)
    
    new_dict = {}
    mstar = subhalo_dict['SubhaloMassInRadType'][:,4].flatten()*(1.0e10)/ilh
    sfids = np.where(mstar >= mstarlim )[0]

    new_dict['SubfindID']=sfids

    for key in subhalo_dict:
        subfield = subhalo_dict[key]
        if key=='count':
            new_dict[key]=subfield
            continue
        new_dict[key]=subhalo_dict[key][sfids]

    #z is mz/mtot

    return new_dict

def make_nlr_bhm():
    o3 = pyfits.open(os.path.expandvars('$HOME/Dropbox/Projects/VELA_sunrise/Sunrise_SPS_Input/Mappings_hard_OIII.fits'))
    o3sed=o3[3].data

    zhard=3 # 3=2x solar
    uhard=3

    #for now, just putting all in brightest of doublets, and only lam < 7000A
    waves=np.asarray([3728.73,4363.15,4068.50,4685.74,4861.32,4958.83,5006.77,5200.00,6087,6300.20,6562.80,6583.34,6716.31])
    inds =np.asarray([0,1,2,3,4,5,6,7,8,9,10,11,12])

    fluxratios=o3sed[2,2,zhard,uhard,inds]
    print(fluxratios)

    bhm_orig = pyfits.open(os.path.expandvars('$HOME/Dropbox/Projects/VELA_sunrise/Sunrise_SPS_Input/newbhmodel_hires.fits'))
    print(bhm_orig.info())
    print(bhm_orig['BHMODEL'].header.cards)
    axes=bhm_orig['AXES']
    print(axes.header.cards)
    sed_orig = bhm_orig['SED'].data
    print(sed_orig.shape)
    print(bhm_orig['SED'].header.cards)
    #sed probably in log W/m ???

    new_sed = np.zeros_like(sed_orig)

    vsigma=50.0 #km/s estimated intrinsic width of lines ??

    lam_m =axes.data['lambda']

    
    loglbol_W=axes.data['log_L_bol']
    for j,llb in enumerate(loglbol_W[0:71]):
        lbol_W=10.0**(llb)
        lo3_W = (u.W)*lbol_W/3500.0 #Heckman+04
        for wl_A,i,fluxrat in zip(waves,inds,fluxratios):
            #need Gaussian model ?
            sigma = (u.m*wl_A/1.0e10)*(vsigma*(u.km/u.s))/astropy.constants.c
            sigma=sigma.to('m')

            lline_W = lo3_W*fluxrat
            
            norm = lline_W*1.0/(((2.0*math.pi)**0.5)*sigma)


            line_func = norm*np.exp(-1.0*((lam_m-wl_A/1.0e10 )**2)/(2.0*sigma.value**2))
            line_sum = np.trapz(line_func,lam_m*u.m)
            
            print(wl_A,lline_W,line_sum, line_sum/lline_W)

            new_sed[:,j]=np.log10(10.0**(new_sed[:,j])  + line_func.value)
            

    bhm_new=copy.copy(bhm_orig)
    bhm_new['SED'].data[:,:]=new_sed

    bhm_new.writeto('bhmodel_nlr_only.fits',overwrite=True)

    return



def plot_simple_bpt(subhalos,radeff=0.1):

    #flux ratios from MAPPINGS.. low U--> SF regions?
    o3 = pyfits.open(os.path.expandvars('$HOME/Dropbox/Projects/VELA_sunrise/Sunrise_SPS_Input/Mappings_hard_OIII.fits'))

    mets=np.asarray([0.25,0.5,1.,2.,4.])

    #just assume solar for AGN component?
    o3sed=o3[3].data
    zhard=3 # 3=2x solar
    uhard=3
    hard_hbo3 = o3sed[2,2,zhard,uhard,4]  #hb=index 4
    hard_hao3 = o3sed[2,2,zhard,uhard,10]  #ha=index 10
    hard_n2o3 = o3sed[2,2,zhard,uhard,11]
    hard_s2o3 = o3sed[2,2,zhard,uhard,12]  #sulfer II line might be useful here?
    
    hard_n2ha = hard_n2o3/hard_hao3
    hard_o3hb = 1.0/hard_hbo3

    #interpolate onto subhalo metallicities
    subz=subhalos['SubhaloGasMetallicitySfrWeighted']
    subzsol=subz/0.02

    #from Kewley+ 2001 n=10, q=2e7, K models
    softmets = np.asarray(['.05','0.2','0.4','1.0','2.0'])
    o3hb=np.zeros_like(softmets,dtype=np.float64)
    hahb=np.zeros_like(softmets,dtype=np.float64)
    n2hb=np.zeros_like(softmets,dtype=np.float64)

    for i,sm in enumerate(softmets):
        o3hb_Z,hahb_Z,n2hb_Z=get_line_fluxes(Z=sm)
        o3hb[i]=o3hb_Z
        hahb[i]=hahb_Z
        n2hb[i]=n2hb_Z

    #o3hb = np.asarray([0.2814,0.9891,0.8876,0.3004,2.3775e-2])
    #hahb = np.asarray([2.905,2.875,2.868,2.945,3.046])
    #n2hb = np.asarray([5.17e-2,0.1666,0.4517,1.220,1.329])

    #usoft = -4
    #soft_hbo3 = np.interp(subzsol,mets,o3sed[0,3,:,usoft,4])  #hb=index 4
    #soft_hao3 = np.interp(subzsol,mets,o3sed[0,3,:,usoft,10])  #ha=index 10
    #soft_n2o3 = np.interp(subzsol,mets,o3sed[0,3,:,usoft,11])
    #soft_s2o3 = np.interp(subzsol,mets,o3sed[0,3,:,usoft,12])  #sulfer II line might be useful here?
    
    #soft_n2ha = soft_n2o3/soft_hao3
    #soft_o3hb = 1.0/soft_hbo3
    soft_o3hb=np.interp(subzsol,np.float64(softmets),o3hb)
    soft_hahb=np.interp(subzsol,np.float64(softmets),hahb)
    soft_n2hb=np.interp(subzsol,np.float64(softmets),n2hb)

    soft_n2ha=soft_n2hb/soft_hahb
    

    bh_lbol= radeff*((u.Msun/u.year)*subhalos['SubhaloBHMdot']*((1.0e10)/ilh)/(0.978*1.0e9/ilh))*(astropy.constants.c**2)
    bh_lbol_lsun = bh_lbol.to('solLum')

    bh_lo3_lsun = bh_lbol_lsun/3500.0 #Heckman+ 04
    bh_lhb_lsun = bh_lo3_lsun*hard_hbo3
    bh_lha_lsun = bh_lo3_lsun*hard_hao3
    bh_ln2_lsun = bh_lo3_lsun*hard_n2o3

    sf_lha_ergs = (u.erg/u.s)*np.float64(subhalos['SubhaloSFRinRad'])/(7.9e-42)  #Kennicutt 98

    sf_lha_lsun = sf_lha_ergs.to('solLum')

    #sf_lo3_lsun = sf_lha_lsun/soft_hao3
    #sf_lhb_lsun = sf_lo3_lsun*soft_hbo3
    #sf_ln2_lsun = sf_lo3_lsun*soft_n2o3
    sf_lhb_lsun = sf_lha_lsun/soft_hahb
    sf_lo3_lsun = sf_lhb_lsun*soft_o3hb
    sf_ln2_lsun = sf_lhb_lsun*soft_n2hb
    
    total_o3 = (sf_lo3_lsun + bh_lo3_lsun)
    total_hb = (sf_lhb_lsun + bh_lhb_lsun)
    total_ha = (sf_lha_lsun + bh_lha_lsun)
    total_n2 = (sf_ln2_lsun + bh_ln2_lsun)

    total_o3hb = total_o3/total_hb
    total_n2ha = total_n2/total_ha


    mstar_msun = subhalos['SubhaloMassInRadType'][:,4]*(1.0e10)/ilh



    #from Hopkins+ 2007... Lbol/Lhard
    #Lbol/Lband=c1*(L/1e10 Lsun)^k1 + c2*(L/1e10 Lsun)^k2
    #for 2-10kev, (10.83,0.28,6.08,-0.02) = (c1,k1,c2,k2)
    c1,k1,c2,k2=10.83,0.28,6.08,-0.02

    Lbol_Lhard=c1*(bh_lbol_lsun.value/1.0e10)**k1 + c2*(bh_lbol_lsun.value/1.0e10)**k2
    Lhard_lsun = bh_lbol_lsun.value/Lbol_Lhard

    #print(np.median(Lhard_lsun),np.max(Lhard_lsun))

    xi = np.where(total_o3hb.value > 2.0)[0]
    xi = np.where(total_o3.value > 10.0**(9.0))[0]

    xi2=np.where(bh_lbol_lsun.value > 1.0e11)[0]

    xi3 = np.where(total_o3.value > (10.0**8.0) )[0]
    
    #print(xi.shape, xi2.shape, xi3.shape, bh_lbol_lsun.shape)

    asi=np.argsort(bh_lbol_lsun.value[xi])
    asi=np.argsort(total_o3.value[xi])
    
    #print(subhalos['SubfindID'][xi[asi]], total_o3hb.value[xi[asi]], np.log10(bh_lbol_lsun.value[xi[asi]]) )


    
    lm = np.log10(mstar_msun)-0.51
    #y = 290.2 − 76.34m + 6.69m2 − 0.1955m3
    div=290.2 - 76.34*lm + 6.69*lm**2 - 0.1955*lm**3
    mxi=np.where(np.log10(total_o3hb) > div)[0]
    #print(mxi.shape)



    print('Mstar: ', total_o3.shape)
    print('OIII:  ', xi3.shape)
    print('Lbol:  ', xi2.shape)
    print('Mex:   ', mxi.shape)
    

    bothi = np.where(np.logical_and(np.log10(total_o3hb) > div, bh_lbol_lsun.value > 1.0e11) )[0]
    #print('both ',bothi.shape)
    
    plot_filen = 'illustris_simple_bpt.pdf'
    f1 = pyplot.figure(figsize=(4.0,3.0), dpi=300)
    pyplot.subplots_adjust(left=0.17, right=0.98, bottom=0.14, top=0.98,wspace=0.0,hspace=0.0)
    axi=f1.add_subplot(1,1,1)


    #this is where the raw AGN emission spectrum sits
    axi.plot(np.log10(hard_n2ha),np.log10(hard_o3hb),'*',markersize=15,color='Orange')


    #softer spectrum estimates
    axi.plot(np.log10(soft_n2ha),np.log10(soft_o3hb),'o',markersize=2,color='Blue')


    #combined spectrum estimates
    axi.plot(np.log10(total_n2ha),np.log10(total_o3hb),'s',markersize=2,color='Black')

    axi.plot(np.log10(total_n2ha[mxi]),np.log10(total_o3hb[mxi]),'s',markersize=3,color='Green')


    sp=10

    axi.set_ylim(-1,1.5)
    axi.set_xlim(-2,0.5)
    axi.set_ylabel(r'$[O III]/H \beta $',size=sp)
    axi.set_xlabel(r'$[N II]/H \alpha $',size=sp)
    axi.tick_params(labelsize=sp)

    f1.savefig(plot_filen,dpi=300)
    pyplot.close(f1)



    
    plot_filen = 'illustris_simple_mex.pdf'
    f1 = pyplot.figure(figsize=(3.5,3.0), dpi=300)
    pyplot.subplots_adjust(left=0.20, right=0.98, bottom=0.14, top=0.98,wspace=0.0,hspace=0.0)
    axi=f1.add_subplot(1,1,1)


    #this is where the raw AGN emission spectrum sits
    #axi.plot(np.log10(1.0e12),np.log10(hard_o3hb),'*',markersize=15,color='Orange')
    axi.plot([9.0,13.0], [np.log10(hard_o3hb),np.log10(hard_o3hb)],marker=None,linestyle='dashed',color='Orange')

    #softer spectrum estimates
    axi.plot(np.log10(mstar_msun),np.log10(soft_o3hb),'o',markersize=1,color='Blue')


    #combined spectrum estimates
    axi.plot(np.log10(mstar_msun),np.log10(total_o3hb),'s',markersize=1,color='Black')
    axi.plot(np.log10(mstar_msun[xi2]),np.log10(total_o3hb[xi2]),'^',markersize=2,color='Red')
    axi.plot(np.log10(mstar_msun[mxi]),np.log10(total_o3hb[mxi]),'o',markersize=2,color='Green')

    axi.annotate('pure AGN',[10.5,1.25],ha='center',va='center',size=10,color='Orange')

    sp=10

    axi.set_ylim(-1,1.5)
    axi.set_xlim(9.8,12.1)
    axi.set_ylabel(r'$[O III]/H \beta $',size=sp)
    axi.set_xlabel(r'$M_*$',size=sp)
    axi.tick_params(labelsize=sp)

    f1.savefig(plot_filen,dpi=300)
    pyplot.close(f1)



    d = illcos.luminosity_distance(1.5)
    
    hb_flx= (total_hb/(4.0*math.pi*d**2) ).to('erg/(s cm^2)')
    o3_flx= (total_o3/(4.0*math.pi*d**2) ).to('erg/(s cm^2)')
    o3_flx_sf = (sf_lo3_lsun/(4.0*math.pi*d**2) ).to('erg/(s cm^2)')

    l3d=np.log10(2.0e-17)
    
    plot_filen = 'illustris_simple_fluxes.pdf'
    f1 = pyplot.figure(figsize=(3.5,3.0), dpi=300)
    pyplot.subplots_adjust(left=0.20, right=0.98, bottom=0.18, top=0.98,wspace=0.0,hspace=0.0)
    axi=f1.add_subplot(1,1,1)
    axi.locator_params(nbins=5,prune='both')
    

    axi.plot(np.log10(hb_flx.value),np.log10(o3_flx_sf.value),'o',markersize=1,color='Blue')


    axi.plot(np.log10(hb_flx.value),np.log10(o3_flx.value),'s',markersize=1,color='Black')
    axi.plot(np.log10(hb_flx.value[xi2]),np.log10(o3_flx.value[xi2]),'^',markersize=2,color='Red')
    axi.plot(np.log10(hb_flx.value[mxi]),np.log10(o3_flx.value[mxi]),'o',markersize=2,color='Green')

    axi.plot([l3d,l3d],[-20,-10],marker=None,linestyle='dotted',color='Gray')

    pyplot.legend(['SF only','SF+AGN','$L_{BH}/L_{\odot} > 10^{11}$','MEx','3D-HST $3\sigma$'],loc='upper left',fontsize=10,ncol=2,markerscale=4,framealpha=1.0,handletextpad=0.5,columnspacing=0.5)

    axi.plot([-20,-10],[l3d,l3d],marker=None,linestyle='dotted',color='Gray')
    
    

    sp=10

    axi.set_ylim(-17.1,-14.1)
    axi.set_xlim(-17.1,-14.9)
    axi.set_xlabel(r'$f(H \beta)$ $[erg/s/cm^2]$',size=sp)
    axi.set_ylabel(r'$f([O III])$ $[erg/s/cm^2]$',size=sp)
    axi.tick_params(labelsize=sp)

    f1.savefig(plot_filen,dpi=300)
    pyplot.close(f1)




    


    return #bh_lbol_lsun, bh_lo3_lsun, sf_lo3_lsun, sf_lha_lsun


if __name__=="__main__":


    subhalos=get_subhalos()
    
