import os
import sys
import numpy as np
import glob
import gfs_sublink_utils as gsu
import shutil
import math
import astropy
import astropy.io.fits as pyfits
#import matplotlib
#import matplotlib.pyplot as plt
import scipy
import scipy.ndimage
#import make_color_image
import numpy.random as random
import congrid
import tarfile
import string
import astropy.io.ascii as ascii



#construct real illustris lightcones from individual images

#in parallel, produce estimated Hydro-ART surveys based on matching algorithms -- high-res?


def build_lightcone_images(image_info_file,run_type='images'):

    data=ascii.read(image_info_file)
    print(data)

    full_npix=data['full_npix'][0]

    #get expected shape
    test_file=os.path.join(data['run_dir'][0],'broadbandz.fits')
    tfo =pyfits.open(test_file)
    print(tfo.info())
    cube=tfo['CAMERA0-BROADBAND-NONSCATTER'].data
    cubeshape=cube.shape
    print(cubeshape)
    auxcube=tfo['CAMERA0-AUX'].data

    filters_hdu = tfo['FILTERS']

    lightcone_dir=os.path.abspath(os.path.dirname(image_info_file))
    print('Constructing lightcone data from: ', lightcone_dir)

    output_dir = os.path.join(lightcone_dir,image_info_file.rstrip('.txt'))
    print('Saving lightcone outputs in: ', output_dir)

    N_filters = cubeshape[0]
    N_aux=auxcube.shape[0]

    image_cube = np.zeros((N_filters,full_npix,full_npix),dtype=np.float64)
    aux_cube = np.zeros((N_aux,full_npix,full_npix),dtype=np.float64)

    for origin_i,origin_j,run_dir,this_npix in zip(data['origin_i'],data['origin_j'],data['run_dir'],data['this_npix']):
        this_cube = pyfits.open(os.path.join(run_dir,'broadbandz.fits'))['CAMERA0-BROADBAND-NONSCATTER'].data
        i_tc=0
        j_tc=0
        i_tc1=this_npix
        j_tc1=this_npix

        if origin_i < 0:
            i0=0
            i_tc=-1*origin_i
        if origin_j < 0:
            j0=0
            j_tc=-1*origin_j
        if i0+this_npix > full_npix:
            i1=full_npix
            i_tc1= full_npix-i0   #this_npix - (i0+this_npix-full_npix)
        if j0+this_npix > full_npix:
            j1=full_npix
            j_tc1= full_npix-j0

        sub_cube1=image_cube[:,i0:i1,j0:j1]
        this_subcube=this_cube[:,i_tc:i_tc1,j_tc:j_tc1]
        print(this_cube.shape, this_npix, sub_cube1.shape, this_subcube.shape)

        image_cube[:,i0:i1,j0:j1] = sub_cube1 + this_subcube



    #convert units before saving.. or save both?

    return
