import numpy as np
import astropy
from astropy.io import fits
from astropy import wcs
import copy

def blank_image(header_only=True,scale_arcsec=0.5,w_deg=0.30,h_deg=0.30,crval1=0.0,crval2=0.0):



    w_np=np.int64(w_deg/(scale_arcsec/3600.0))
    h_np=np.int64(h_deg/(scale_arcsec/3600.0))

    blankwcs=wcs.WCS(naxis=2)
    blankwcs.wcs.crpix=[int(w_np)/2,int(h_np)/2]
    blankwcs.wcs.cdelt=[scale_arcsec/3600,scale_arcsec/3600]
    blankwcs.wcs.crval=[ crval1,   crval2]
    blankwcs.wcs.ctype=['RA---TAN','DEC--TAN']


    basic_cd=np.array([[ scale_arcsec/3600,  0.0],[ 0.0,   scale_arcsec/3600]])
    blankwcs.wcs.cd=basic_cd
    blank_wcs=copy.copy(blankwcs)
    blank_header=blank_wcs.to_header()

    blank_header['PIXSCALE']=(scale_arcsec,'arcsec')
    blank_header['W_DEG']=(w_deg,'degrees')
    blank_header['H_DEG']=(h_deg,'degrees')

    if header_only is False:
        blank_image=np.ndarray(shape=(h_np,w_np),dtype=np.float32)
        imhdu=fits.PrimaryHDU(blank_image,header=blank_header)
    else:
        imhdu=fits.PrimaryHDU(header=blank_header)


    return imhdu
