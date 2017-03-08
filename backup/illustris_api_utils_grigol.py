import cProfile
import pstats
import math
import string
import sys
import struct
import numpy as np
import cPickle
import asciitable
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import scipy.odr as odr
import glob
import os
import make_color_image
import make_fake_wht
import gzip
import tarfile
import shutil
import cosmocalc
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
from multiprocessing import Process, Queue, current_process
import time
import illustris_lightcone_catalogs as ilc
import translate_coordinates as tc

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

    #subhalo metadata ?
    sub = get(url)
    #print sub['pos_x']
    sub_pos_x = sub['pos_x']
    sub_pos_y = sub['pos_y']
    sub_pos_z = sub['pos_z']

    if os.path.lexists(checkf) and clobber is False:
        if verbose:
            print "Subhalo cutout exists, skipping: ", checkf
        return checkf, (sub_pos_x,sub_pos_y,sub_pos_z)
    else:
        if verbose:
            print "Getting subhalo cutout: ", checkf
        file = get(url+"/cutout.hdf5",params,savepath)
        return file, (sub_pos_x,sub_pos_y,sub_pos_z)


def process_subhalo(sim,snap,sfid,camera_obj,params=defaultparams,savepath=None,verbose=True,clobber=False,getlabel='StellarMass',resample=False):
    #defaultparams={'stars':'Coordinates,Velocities,GFM_StellarFormationTime,GFM_Metallicity,GFM_InitialMass,Masses'}

    file, pos = get_subhalo(sim,snap,sfid,params=defaultparams,savepath=savepath,verbose=verbose,clobber=clobber)
    substuff = None

    #load HDF5 file relevant fields
    with h5py.File(file) as h5f:
        masses = np.asarray(h5f['PartType4']['Masses'])*(1.0e10)/ilh  #convert to solar masses
        formfactors = np.asarray(h5f['PartType4']['GFM_StellarFormationTime'])  #in scale factor
        xpos = np.asarray(h5f['PartType4']['Coordinates'][:,0])-pos[0] #in ckpc/h
        ypos = np.asarray(h5f['PartType4']['Coordinates'][:,1])-pos[1]
        zpos = np.asarray(h5f['PartType4']['Coordinates'][:,2])-pos[2]

    print masses.shape, ypos.shape, np.mean(ypos)


    #(optional) resample star particles
    
    #transform into lightcone camera coords

    #project into an image to return
    #ideas for options:  F435W, F606W, F850LP, F105W, F125W, F160W, NIRCAM, Mstar, Mgas, SFR, Simonsized H-alpha lines?

    return (file, substuff)


def worker(input,output,**kwargs):
    for func, args in iter(input.get,'STOP'):
        f = calculate(func,args,**kwargs)
        output.put(f)

def calculate(func,args,**kwargs):
    result = func(*args,**kwargs)
    #return '%s was given %s and got %s' % \
    #         (current_process().name, args, result)
    return result

def get_lightcone_images_threaded(lcfile,geofile,sim='Illustris-2',clobber=False,savepath='/astro/snyder_lab2/Illustris',Np=2,maxq=10000,lim=None):
    data = ascii.read(lcfile)
    snapnums = np.asarray(data['col1'],dtype='str')
    sfids = np.asarray(data['col2'],dtype='str')
    x_cmpc = np.asarray(data['col13'])
    y_cmpc = np.asarray(data['col14'])
    z_cmpc = np.asarray(data['col15'])
    coords_ra = np.asarray(data['col3']) #degrees
    coords_dec = np.asarray(data['col4']) #degrees

    catalog = ilc.process_lightcone_catalog(lightcone=geofile,basedir='.')
    print catalog.delb_arcmin


    NUMBER_OF_PROCESSES=Np
    task_queue = Queue()
    done_queue = Queue()
    TASKS = []
    TASKS_DONE = []
    TASKS_LEFT = []

    N_objects = snapnums.shape[0]
    #figure out how to drain the queue before it fills the pipe

    if lim is None:
        lim=N_objects

    for i,sn in enumerate(snapnums[0:lim]):
        this_sfid = sfids[i]
        this_x = x_cmpc[i]
        this_y = y_cmpc[i]
        this_z = z_cmpc[i]
        this_ra = coords_ra[i]
        this_dec = coords_dec[i]

        task = (process_subhalo,(sim,sn,this_sfid,None))
        if i <= maxq:
            task_queue.put(task)
            TASKS.append(task)
        else:
            TASKS_LEFT.append(task)

    for p in range(NUMBER_OF_PROCESSES):
        Process(target=worker,args=(task_queue,done_queue),kwargs={'savepath':savepath,'verbose':True,'clobber':clobber}).start()

    cutout_files = []        

    while len(TASKS_LEFT) > 0:
        cutout_files.append(done_queue.get())
        newtask = TASKS_LEFT.pop()
        task_queue.put(newtask)

    for i in range(min(maxq,lim)):
        cutout_files.append(done_queue.get())


    #build up images/analysis here!
    print cutout_files[0:5]
    print cutout_files[-5:]
    print len(cutout_files)

    for p in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')

    return None





def get_lightcone_images(lcfile,geofile,sim='Illustris-2',clobber=False,savepath='/astro/snyder_lab2/Illustris'):

    data = ascii.read(lcfile)
    snapnums = np.asarray(data['col1'],dtype='str')
    sfids = np.asarray(data['col2'],dtype='str')


    for i,sn in enumerate(snapnums):
        this_sfid = sfids[i]
        #get file.  will skip download if it already exists
        f = get_subhalo(sim,sn,this_sfid,savepath=savepath,verbose=True,clobber=clobber)
        #obtain fields, place at desired position, project, compute densities and luminosities
        
        #for projection, how?  use lightcone direction? if so, must save somewhere!

    return


def do_lightcone_images():

    lcfile = '/astro/snyder_lab2/Illustris/Lightcones/Illustris-2_RADEC_hudfwide_75Mpc_7_6_xyz_corners.txt'
    geofile = '/astro/snyder_lab2/Illustris/Lightcones/hudfwide_75Mpc_7_6_fixedh_xyz_NEW.txt'
    savepath = '/astro/snyder_lab2/Illustris'

    #note, for Gordon, stage lightcone data on local scratch!!

    #st1 = time.time()
    #get_lightcone_images(lcfile,geofile,sim='Illustris-2',clobber=False,savepath=savepath)
    #et1 = time.time()

    st = time.time()
    result = get_lightcone_images_threaded(lcfile,geofile,sim='Illustris-2',clobber=True,savepath='/astro/snyder_lab2/Illustris',Np=4,lim=25)
    et = time.time()

    print 'Threaded calculation took: ', et-st, ' seconds'
    #print 'Serial calculation took: ', et1-st1, ' seconds'


