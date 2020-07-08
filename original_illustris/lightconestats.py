import math
import string
import sys
import struct
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.colors as pycolors
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import asciitable
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import scipy.odr as odr
from scipy.integrate import quad
import glob
import os
import gzip
import shutil
import astropy.io.ascii as ascii
import warnings
import subprocess
import astropy
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.cosmology import WMAP7,z_at_value
import copy
import datetime
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import PolyCollection
import PyML
import PyML.machinelearning
import gfs_twod_histograms as gth
import pandas
import h5py
import make_color_image
import photutils




def plot_mag_versus_z():



    lcfiles=np.asarray(glob.glob('/astro/snyder_lab2/Illustris/Lightcones/CEERS/Illustris-1_RADEC_ceers_75Mpc_9_8_???.txt'))


    fig=pyplot.figure(figsize=(6,3)) 
    fig.subplots_adjust(bottom=0.15,left=0.09)
    axi=fig.add_subplot(1,1,1)

    for lcf in lcfiles:
        dat=ascii.read(lcf)
        mag=dat['col39']
        z=dat['col9']
        print(mag.shape)
        print(z.shape)
        axi.plot(z,mag,'o')
    
    axi.annotate('3x28 arcmin$^2$',(0.75,0.10),xycoords='axes fraction',color='black',ha='center',va='center',size=9)

    axi.set_xlim(6.0,12.0)
    axi.set_ylim(26.0,35.0)
    axi.set_xlabel('redshift')
    axi.set_ylabel('rest-g apparent AB mag')

    fig.savefig('Illustris-1_9_8_mag_vs_z.pdf')
    pyplot.close(fig)




    lcfiles=np.asarray(glob.glob('/astro/snyder_lab2/Illustris/Lightcones/CEERS/Illustris-1_RADEC_ceers_75Mpc_11_10_???.txt'))
    print(lcfiles.shape)

    fig=pyplot.figure(figsize=(6,3)) 
    fig.subplots_adjust(bottom=0.15,left=0.09)
    axi=fig.add_subplot(1,1,1)

    for lcf in lcfiles:
        dat=ascii.read(lcf)
        mag=dat['col39']
        z=dat['col9']
        print(mag.shape)
        print(z.shape)
        axi.plot(z,mag,'o')
    
    axi.annotate('3x8 arcmin$^2$',(0.75,0.10),xycoords='axes fraction',color='black',ha='center',va='center',size=9)

    axi.set_xlim(6.0,12.0)
    axi.set_ylim(26.0,35.0)
    axi.set_xlabel('redshift')
    axi.set_ylabel('rest-g apparent AB mag')

    fig.savefig('Illustris-1_11_10_mag_vs_z.pdf')
    pyplot.close(fig)






    lcfiles=np.asarray(glob.glob('/astro/snyder_lab2/Illustris/Lightcones/CEERS/Illustris-1_RADEC_ceers_75Mpc_7_6_???.txt'))
    print(lcfiles.shape)

    fig=pyplot.figure(figsize=(6,3)) 
    fig.subplots_adjust(bottom=0.15,left=0.09)
    axi=fig.add_subplot(1,1,1)

    for lcf in lcfiles:
        dat=ascii.read(lcf)
        mag=dat['col39']
        z=dat['col9']
        print(mag.shape)
        print(z.shape)
        axi.plot(z,mag,'o')
    
    axi.annotate('3x140 arcmin$^2$',(0.75,0.10),xycoords='axes fraction',color='black',ha='center',va='center',size=9)

    axi.set_xlim(6.0,12.0)
    axi.set_ylim(26.0,35.0)
    axi.set_xlabel('redshift')
    axi.set_ylabel('rest-g apparent AB mag')

    fig.savefig('Illustris-1_7_6_mag_vs_z.pdf')
    pyplot.close(fig)




    return



if __name__=="__main__":
    plot_mag_versus_z()
