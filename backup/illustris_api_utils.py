import cProfile
import pstats
import math
import string
import sys
import struct
import numpy as np
import asciitable
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import scipy.odr as odr
import glob
import os

import gzip
import tarfile
import shutil
import congrid
import astropy.io.ascii as ascii
import warnings
import subprocess
import photutils
import astropy
import astropy.cosmology
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.cosmology import WMAP7,z_at_value
from astropy.coordinates import SkyCoord
import copy
import medianstats_bootstrap as msbs
import illustris_python as ilpy
import h5py
import requests

ilh = 0.704
illcos = astropy.cosmology.FlatLambdaCDM(H0=70.4,Om0=0.2726,Ob0=0.0456)

defaultparams={'stars':'Coordinates,Velocities,GFM_StellarFormationTime,GFM_Metallicity,GFM_InitialMass,Masses'}

baseUrl = 'http://www.illustris-project.org/api/'
headers = {"api-key":"117782db3bf216d7ce7a04d0c9034601"}

def get(path, params=None,savepath=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)
    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()
    
    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically
    
    if 'content-disposition' in r.headers:
        file_basename = r.headers['content-disposition'].split("filename=")[1]

        if savepath is not None:
            filename = os.path.join(savepath,file_basename)
        else:
            filename = file_basename

        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r



def get_subhalo(sim,snap,sfid,params=defaultparams,savepath=None,verbose=True,clobber=False):


    relative_path = sim+"/snapshots/" + str(snap) + "/subhalos/" + str(sfid)

    url = "http://www.illustris-project.org/api/"+relative_path
    #sub = get(url)

    if savepath is not None:
        savepath = os.path.join(savepath,relative_path)
        if not os.path.lexists(savepath):
            os.makedirs(savepath)
    
        checkf = os.path.join(savepath,'cutout_'+str(sfid)+'.hdf5')
    else:
        checkf = 'cutout_'+str(sfid)+'.hdf5'


    if os.path.lexists(checkf) and clobber is False:
        if verbose:
            print("Subhalo cutout exists, skipping: ", checkf)
        return checkf
    else:
        if verbose:
            print("Getting subhalo cutout: ", checkf)
        file = get(url+"/cutout.hdf5",params,savepath)
        return file



def get_lightcone_images(lcfile,sim='Illustris-2'):

    data = ascii.read(lcfile)
    snapnums = np.asarray(data['col1'],dtype='str')
    sfids = np.asarray(data['col2'],dtype='str')

    savepath = '/astro/snyder_lab2/Illustris'

    for i,sn in enumerate(snapnums):
        this_sfid = sfids[i]
        #get file.  will skip download if it already exists
        f = get_subhalo(sim,sn,this_sfid,savepath=savepath,verbose=True,clobber=False)
        #obtain fields, place at desired position, project, compute densities and luminosities
        
        #for projection, how?  use lightcone direction? if so, must save somewhere!

    return


def do_lightcone_images():

    get_lightcone_images('/astro/snyder_lab2/Illustris/Lightcones/Illustris-2_RADEC_hudfwide_75Mpc_7_6_xyz.txt')
