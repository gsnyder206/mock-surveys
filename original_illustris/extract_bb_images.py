import astropy
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
import numpy as np
import make_color_image
import matplotlib
import matplotlib.pyplot as pyplot


def quickrgb(image_file,alph=1.0,Q=1.0,mid=500,delt=500,fov_comoving=True,rgb_factors=[2.0,1.5,1.15],**kwargs):

    ims,fils,hduheader,bbheader=extract_bb_images(**kwargs)

    if fov_comoving is True:
        redshift=bbheader['REDSHIFT']
        zdelt=int(float(delt)/(1.0 + redshift))
    else:
        zdelt=delt
            
    bb=rgb_factors[2]*ims[0][mid-zdelt:mid+zdelt,mid-zdelt:mid+zdelt]
    gg=rgb_factors[1]*ims[1][mid-zdelt:mid+zdelt,mid-zdelt:mid+zdelt]
    rr=rgb_factors[0]*ims[2][mid-zdelt:mid+zdelt,mid-zdelt:mid+zdelt]
    
    
    fig = pyplot.figure(figsize=(6,6), dpi=600)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)
    axi = fig.add_subplot(1,1,1)
    axi.set_xticks([]) ; axi.set_yticks([])
    
    rgbthing = make_color_image.make_interactive_nasa(bb,gg,rr,alph,Q)
    axi.imshow(rgbthing,interpolation='nearest',aspect='auto',origin='lower')
    
    fig.savefig(image_file,dpi=300,facecolor='Black')
    pyplot.close(fig)

    return rgbthing


def extract_bb_images(filename='broadbandz.fits',hdukey='CAMERA0-BROADBAND-NONSCATTER',fils=['jwst/nircam_f115w','jwst/nircam_f150w','jwst/nircam_f200w']):

    bb=pyfits.open(filename)
    hduheader=bb[hdukey].header
    hdudata=bb[hdukey].data

    ims=[]
    
    fils_from_hdu=bb['FILTERS'].data['filter']
    bbhduheader=bb['BROADBAND'].header
    
    for f in fils:
        if f in fils_from_hdu:
            fi=np.argwhere(fils_from_hdu==f)[0][0]
            im=hdudata[fi,:,:]
            ims.append(im)
            
    return ims,fils,hduheader,bbhduheader



def extract_integrated_spectrum(filename,camera=0,ltype='out'):

    bb=pyfits.open(filename)
    iqdata=bb['INTEGRATED_QUANTITIES'].data
    iqhead=bb['INTEGRATED_QUANTITIES'].header

    lambda_m=iqdata['lambda']
    l_Wm=iqdata['L_lambda_'+ltype+str(camera)]
    l_Wm_ns=iqdata['L_lambda_nonscatter'+str(camera)]

    
    return l_Wm, lambda_m, l_Wm_ns

