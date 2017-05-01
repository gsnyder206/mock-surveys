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
from astropy.convolution import *
import copy


if __name__=="__main__":
    files = np.asarray(glob.glob('*.fits'))

    for fn in files:
        
