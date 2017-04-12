import os
import sys
import numpy as np
import glob
import gfs_sublink_utils as gsu
import shutil
import math
import astropy
import astropy.io.fits as pyfits
import matplotlib
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage
import make_color_image
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


    lightcone_dir=os.path.dirname(image_info_file)
    print('Constructing lightcone data from: ', lightcone_dir)

    output_dir = os.path.join(lightcone_dir,image_info_file.rstrip('.txt'))
    print('Saving lightcone outputs in: ', lightcone_dir)

    return
