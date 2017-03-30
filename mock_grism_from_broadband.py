import astropy
import astropy.io.fits as pyfits
import matplotlib
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage
import showgalaxy
import make_color_image
import numpy.random as random
import gfs_sublink_utils as gsu
import numpy as np


def make_simple_trace(bbfile='grism.fits',outname='grismtrace',ybox=None,xbox=None,noisemaxfact=0.05,alph=1.0,Q=1.0):

    go=pyfits.open(bbfile)
    
    
    redshift=go['BROADBAND'].header['REDSHIFT']
    wfc3_pix_as=0.13
    g141_nm_per_pix=4.65
    
    min_lam=1.075
    max_lam=1.700
    
    hdu=go['CAMERA0-BROADBAND-NONSCATTER']
    cube=hdu.data #L_lambda units! 
    #cube=np.flipud(cube) ; print(cube.shape)
    
    fil=go['FILTERS']
    lamb=fil.data['lambda_eff']*1.0e6
    flux=fil.data['L_lambda_eff_nonscatter0']
    
    g141_i = (lamb >= min_lam) & (lamb <= max_lam)
    
    arcsec_per_kpc= gsu.illcos.arcsec_per_kpc_proper(redshift)
    kpc_per_arcsec=1.0/arcsec_per_kpc.value
    
    im_kpc=hdu.header['CD1_1']
    print('pix size kpc: ', im_kpc)
    
    wfc3_kpc_per_pix=wfc3_pix_as*kpc_per_arcsec
    total_width_pix=(1.0e3)*(max_lam-min_lam)/g141_nm_per_pix
    total_width_kpc=total_width_pix*wfc3_kpc_per_pix
    
    total_width_impix=int(total_width_kpc/im_kpc)
    
    delta_lam=(max_lam-min_lam)/total_width_impix  #microns/pix
    
    psf_arcsec=0.18
    psf_kpc=psf_arcsec*kpc_per_arcsec
    psf_impix=psf_kpc/im_kpc
    
    
    imw_cross=200
    imw_disp=total_width_impix+imw_cross
    Np=cube.shape[-1]
    mid = np.int64(Np/2)
    delt=np.int64(imw_cross/2)
    output_image=np.zeros_like( np.ndarray(shape=(imw_disp,imw_cross),dtype='float' ))
    #r = r[mid-delt:mid+delt,mid-delt:mid+delt]
    output_image.shape
    small_cube=cube[g141_i,mid-delt:mid+delt,mid-delt:mid+delt]
    
    for i,l in enumerate(lamb[g141_i]):
        di=int( (l-min_lam)/delta_lam )
        this_cube=small_cube[i,:,:]*l**2  #convert to Janskies-like
        #if i==17:
        #    this_cube[30,30] = 1.0e3
        #print(i,l/(1.0+redshift),int(di),np.sum(this_cube),this_cube.shape,output_image.shape,output_image[di:di+imw_cross,:].shape)
        output_image[di:di+imw_cross,:]=output_image[di:di+imw_cross,:]+this_cube
        
        
    output_image=scipy.ndimage.gaussian_filter(output_image,sigma=[4,psf_impix/2.355])
    
    new_thing = np.transpose(np.flipud(output_image))
    
    
    nr = noisemaxfact*np.max(new_thing)*random.randn(new_thing.shape[0],new_thing.shape[1])
    
    thing=make_color_image.make_interactive(new_thing+nr,new_thing+nr,new_thing+nr,alph=alph,Q=Q)
    thing=1.0-np.fliplr(np.transpose(thing,axes=[1,0,2]))
    

    f=plt.figure(figsize=(25,6))
    axi=f.add_subplot(1,1,1)
    axi.imshow( (thing),aspect='auto',origin='left',interpolation='nearest')
    f.savefig(outname+'.png',dpi=500)
    plt.close(f)

    #[ybox[0]:ybox[1],xbox[0]:xbox[1]]
    #[50:125,120:820,:]

    new_hdu=pyfits.PrimaryHDU(new_thing)
    new_list=pyfits.HDUList([new_hdu])
    new_list.writeto(outname+'.fits',clobber=True)


    return thing, new_thing


