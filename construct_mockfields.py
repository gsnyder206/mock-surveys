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
import np.random as random
import congrid
import tarfile
import string
import astropy.io.ascii as ascii



#construct real illustris lightcones from individual images

#in parallel, produce estimated Hydro-ART surveys based on matching algorithms -- high-res?


def build_lightcone_images(image_info_file):

    data=ascii.read(image_info_file)
    print(data)

    return
