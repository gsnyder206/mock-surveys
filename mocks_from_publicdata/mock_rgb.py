import astropy
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import numpy as np
import matplotlib
import matplotlib.pyplot as pyplot
import make_color_image
import astropy.units as u
from astropy.cosmology import WMAP7,z_at_value
import glob
import sys
import os



def ceers_field_rgb(alph=1.0,Q=1.0,root='ceers_sam_egs-0_',filtemp='NIRCam_F{}W.fits',remake=False,tempname='sam',subset=None):

    savefile=root+'rgb.pdf'
    print(savefile)

    f=pyplot.figure(figsize=(24.0,4.0),dpi=1200)
    pyplot.subplots_adjust(wspace=0.0,hspace=0.0,top=1.0,right=1.0,left=0.00,bottom=0.00)

    ax=f.add_subplot(1,1,1)
    ax.set_xticks([]) ; ax.set_yticks([])



    if not os.path.lexists('temprgb_'+tempname+'.fits') or remake==True:

        #files in nJ
        print('loading images..')
        f115=1.0*fits.open(root+filtemp.format('115'))[0].data  #these are [short, long]
        print('F115W')
        f150=1.0*fits.open(root+filtemp.format('150'))[0].data
        f200=1.0*fits.open(root+filtemp.format('200'))[0].data
        print('F200W')
        f277=f200 #1.0*fits.open(root+filtemp.format('277'))[0].data
        f356=1.0*fits.open(root+filtemp.format('356'))[0].data
        f444=1.0*fits.open(root+filtemp.format('444'))[0].data
        print('F444W', f444.shape)

        b=1.25*(f115)
        g=0.50*(f150+f200)
        r=0.33*(f277+f356+f444)
        print('added')

        print('alph, Q: ', alph, Q)
        if subset is not None:
            rgbthing = make_color_image.make_interactive_nasa(b[0:subset,0:6*subset],g[0:subset,0:6*subset],r[0:subset,0:6*subset],alph,Q)
        else:
            rgbthing = make_color_image.make_interactive_nasa(b,g,r,alph,Q)
        hdu=fits.PrimaryHDU(rgbthing)
        hdu.writeto('temprgb_'+tempname+'.fits',overwrite=True)
        print('made RGB image', rgbthing.shape)
    else:
        print('loading temprgb_'+tempname+'.fits')
        rgbthing=fits.open('temprgb_'+tempname+'.fits')[0].data
        print(rgbthing.shape)

    #expects [short, long, 3] to match figure dims
    thing=np.transpose(rgbthing,axes=(1,0,2))
    print(thing.shape)

    ax.imshow(thing,interpolation='nearest',aspect='auto',origin='lower')
    print('imshowed')

    f.savefig(savefile,dpi=1200)

    return


if __name__=="__main__":
    if len(sys.argv)>1:
        alph=float(sys.argv[1])
        Q=float(sys.argv[2])
        remake=bool(int(sys.argv[3]))
    else:
        alph=0.3
        Q=4.0
        remake=True

    ceers_field_rgb(alph=alph,Q=Q,remake=remake,root='ceers_sam_egs-2_',tempname='sam',subset=None)
    ceers_field_rgb(alph=alph,Q=Q,remake=remake,root='ceers_hydro_egs-2_',filtemp='jwst_f{}w.fits',tempname='hydro',subset=None)
