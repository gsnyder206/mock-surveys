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

flh = 0.695
fogcos = astropy.cosmology.FlatLambdaCDM(H0=70.4,Om0=0.285,Ob0=0.2850-0.2389)

def foggie_time_series(label,alph=1e3,Q=15,label_redshift=True,label_time=True,psf_fwhm_arcsec=None,fov_comoving=True):
    savefile='timeseries'+label+'.pdf'
    
    f=pyplot.figure(figsize=(12.0,8.0),dpi=500) 
    pyplot.subplots_adjust(wspace=0.0,hspace=0.0,top=1.0,right=1.0,left=0.00,bottom=0.00)
    


    flist=np.sort(np.asarray(glob.glob('RD00??/*_sunrise/input/images/broadbandz.fits')))
    print(flist.shape,flist)
    
    mid=500
    delt=500

    Nx=3
    Ny=2
    Ntot=Nx*Ny
    
    for i,fn in enumerate(flist[0:Ntot]):
        ims,fils,header,bbheader=ebi.extract_bb_images(fn,
                                                       fils=['hst/acs_f606w','jwst/nircam_f115w','jwst/nircam_f277w'],
                                                       hdukey='CAMERA5-BROADBAND-NONSCATTER')

        ax=f.add_subplot(Ny,Nx,i+1)
        ax.set_xticks([]) ; ax.set_yticks([])
        print(ims[2].shape)

        if fov_comoving is True:
            redshift=bbheader['REDSHIFT']
            zdelt=int(float(delt)/(1.0 + redshift))
        else:
            zdelt=delt
            
        b=(0.61**1.0)*ims[0][mid-zdelt:mid+zdelt,mid-zdelt:mid+zdelt]
        g=(1.15**1.0)*ims[1][mid-zdelt:mid+zdelt,mid-zdelt:mid+zdelt]
        r=(2.77**1.0)*ims[2][mid-zdelt:mid+zdelt,mid-zdelt:mid+zdelt]
    


        if psf_fwhm_arcsec is not None:
            redshift=bbheader['REDSHIFT']
            kpc_per_arcsec=fogcos.kpc_proper_per_arcmin(redshift).value/60.0
            kpc_per_pixel=header['CD1_1']
            psf_fwhm_pixels=psf_fwhm_arcsec*kpc_per_arcsec/kpc_per_pixel
            fwhm_pixels=np.asarray([psf_fwhm_pixels,psf_fwhm_pixels,psf_fwhm_pixels])
        else:
            fwhm_pixels=np.asarray([0.0,0.0,0.0])


        alph=alph ; Q=Q
        rgbthing = make_color_image.make_interactive(b,g,r,alph,Q,fwhm_pixels=fwhm_pixels)
        ax.imshow(rgbthing,interpolation='nearest',aspect='auto',origin='lower')
        
        
        if label_redshift is True:
            ax.annotate('z={:6.2f}'.format(bbheader['REDSHIFT']),
                        (0.50,0.85),ha='center',xycoords='axes fraction',
                        color='white',fontsize=10,va='center')
        if label_time is True:
            ax.annotate('t={:6.2f} Gyr'.format(fogcos.age(bbheader['REDSHIFT']).value ),
                        (0.50,0.95),ha='center',xycoords='axes fraction',
                        color='white',fontsize=10,va='center')
        
        #ax.annotate(label[i],(0.50,0.85),ha='center',xycoords='axes fraction',color='white',fontsize=10,va='center',bbox=dict(boxstyle='round', fc="gray", ec="blue"))
        
    f.savefig(savefile,dpi=500)

if __name__=="__main__":
    if len(sys.argv)==2:
        label=sys.argv[1]
        psf_arcsec=None
    elif len(sys.argv)==3:
        label=sys.argv[1]
        psf_arcsec=float(sys.argv[2])
    else:
        label=''
        psf_arcsec=None
        
    foggie_time_series(label,psf_fwhm_arcsec=psf_arcsec)
