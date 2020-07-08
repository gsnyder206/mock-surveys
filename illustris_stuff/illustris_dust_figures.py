import astropy
import astropy.io.ascii as ascii
import numpy as np
import extract_bb_images as ebi
import matplotlib
import matplotlib.pyplot as pyplot
import make_color_image
import astropy.units as u
from astropy.cosmology import WMAP7,z_at_value
import glob
import sys
import os


ilh = 0.704
illcos = astropy.cosmology.FlatLambdaCDM(H0=70.4,Om0=0.2726,Ob0=0.0456)

def plot_image(ax,ims,alph=1e3,Q=15,mid=400,delt=400):
    ax.set_xticks([]) ; ax.set_yticks([])

    b=(0.90**1.0)*ims[0][mid-delt:mid+delt,mid-delt:mid+delt]
    g=(0.45**1.0)*ims[1][mid-delt:mid+delt,mid-delt:mid+delt]
    r=(2.00**1.0)*ims[2][mid-delt:mid+delt,mid-delt:mid+delt]
    
    alph=alph ; Q=Q
    rgbthing = make_color_image.make_interactive(b,g,r,alph,Q)
    ax.imshow(rgbthing,interpolation='nearest',aspect='auto',origin='lower')
        
    return


def illustris_dust_plots(alph=1e3,Q=20):

    source_dir='/astro/snyder_lab2/Illustris/HST26'
    figures_dir='/Users/gsnyder/Dropbox/HST_Cycle26/HST26_Vogelsberger_Dust'

    bb_mw=os.path.join(source_dir,'sh13290/images/broadbandz.fits')
    bb_1=os.path.join(source_dir,'sh13290/images1/broadbandz.fits')
    bb_2=os.path.join(source_dir,'sh13290/images2/broadbandz.fits')

    
    f=pyplot.figure(figsize=(8.0,8.0),dpi=500) 
    pyplot.subplots_adjust(wspace=0.05,hspace=0.05,top=1.0,right=1.0,left=0.1,bottom=0.1)
    
    
    mid=400
    delt=200


    ims_mw,fils_mw,header_mw,bbheader_mw=ebi.extract_bb_images(bb_mw,
                                                               fils=['hst/acs_f435w','hst/acs_f814w','jwst/nircam_f200w'],
                                                               hdukey='CAMERA2-BROADBAND')
    ims_1,fils_1,header_1,bbheader_1=ebi.extract_bb_images(bb_1,
                                                           fils=['hst/acs_f435w','hst/acs_f814w','jwst/nircam_f200w'],
                                                           hdukey='CAMERA2-BROADBAND')
    ims_2,fils_2,header_2,bbheader_2=ebi.extract_bb_images(bb_2,
                                                           fils=['hst/acs_f435w','hst/acs_f814w','jwst/nircam_f200w'],
                                                           hdukey='CAMERA2-BROADBAND')
    
    
    ax= pyplot.subplot2grid((2, 3), (0, 0))
    plot_image(ax,ims_mw,alph=alph,Q=Q,mid=mid,delt=delt)
    ax.annotate('MW dust',(0.5,0.1),xycoords='axes fraction',size=18,color='white',ha='center',va='center')
    
    ax= pyplot.subplot2grid((2, 3), (0, 1))
    plot_image(ax,ims_1,alph=alph,Q=Q,mid=mid,delt=delt)
    ax.annotate('Alt dust 1',(0.5,0.1),xycoords='axes fraction',size=18,color='white',ha='center',va='center')

    ax= pyplot.subplot2grid((2, 3), (0, 2))
    plot_image(ax,ims_2,alph=alph,Q=Q,mid=mid,delt=delt)
    ax.annotate('Alt dust 2',(0.5,0.1),xycoords='axes fraction',size=18,color='white',ha='center',va='center')

    

    l_Wm_mw,lambda_m_mw,l_Wm_mw_ns=ebi.extract_integrated_spectrum(bb_mw,camera=0,ltype='scatter')
    l_Wm_1,lambda_m_1,l_Wm_1_ns=ebi.extract_integrated_spectrum(bb_1,camera=0,ltype='scatter')
    l_Wm_2,lambda_m_2,l_Wm_2_ns=ebi.extract_integrated_spectrum(bb_2,camera=0,ltype='scatter')
   
    ax3=pyplot.subplot2grid((2, 3), (1,0), colspan=3)

    mw,=ax3.loglog(lambda_m_mw*1.0e6,l_Wm_mw,color='black',linestyle='solid',lw=3)
    ns,=ax3.loglog(lambda_m_mw*1.0e6,l_Wm_mw_ns,color='gray',linestyle='solid',lw=1)

    alt1,=ax3.loglog(lambda_m_1*1.0e6,l_Wm_1,color='purple',linestyle='dashed',lw=3)
    alt2,=ax3.loglog(lambda_m_2*1.0e6,l_Wm_2,color='red',linestyle='dotted',lw=3)

    #ax3.loglog(lambda_m_smc*1.0e6,l_Wm_smc_ns,color='pink',linestyle='dashed',lw=1)
    ax3.legend((ns,mw,alt1,alt2),('No dust', 'MW', 'Alt1', 'Alt2'),loc='upper right', fontsize=18)
    
    ls=18
    
    ax3.set_xlim(0.1,2.0)
    ax3.set_ylim(1.0e42,3.0e45)
    ax3.set_xlabel('wavelength (microns)',size=ls)
    ax3.set_ylabel(r'$F_{\lambda}$',size=ls)
    
    savefile=os.path.join(figures_dir,'illustris_sn068_sh13290_dustmodels.pdf')
    f.savefig(savefile,dpi=500)

    

if __name__=="__main__":


    illustris_dust_plots()
