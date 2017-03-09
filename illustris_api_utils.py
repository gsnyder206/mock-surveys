import cProfile
import pstats
import math
import string
import sys
import struct
import numpy as np
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import scipy.odr as odr
import glob
import os
import make_color_image
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
#import illustris_python as ilpy
import h5py
import requests
from multiprocessing import Process, Queue, current_process
import time
#import illustris_lightcone_catalogs as ilc
#import translate_coordinates as tc
import gfs_sublink_utils as gsu


ilh = 0.704
illcos = astropy.cosmology.FlatLambdaCDM(H0=70.4,Om0=0.2726,Ob0=0.0456)

start_time = time.time()

defaultparams={'stars':'ParticleIDs','Coordinates,Velocities,GFM_StellarFormationTime,GFM_Metallicity,GFM_InitialMass,Masses','gas':'Coordinates,Density,ElectronAbundance,Masses,Velocities,Volume,SubfindDensity,Potential,InternalEnergy,StarFormationRate,GFM_Metallicity,GFM_AGNRadiation,GFM_WindDMVelDisp,GFM_CoolingRate,NeutralHydrogenAbundance,SmoothingLength,SubfindHsml,SubfindVelDisp,NumTracers,ParticleIDs','dm':'Coordinates,Velocities,Potential','bhs':'all'}


baseUrl = 'http://www.illustris-project.org/api/'
headers = {"api-key":"117782db3bf216d7ce7a04d0c9034601"}


class subhalo_observation:
    def __init__(self,sim,sn,sfid,i,label):
        self.sim=sim
        self.sn=sn
        self.sfid=sfid
        self.i=i
        self.label=label



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

    sim_url="http://www.illustris-project.org/api/"+sim

    snap_url=sim_url+"/snapshots/"+str(snap)


    if savepath is not None:
        savepath = os.path.join(savepath,relative_path)
        if not os.path.lexists(savepath):
            os.makedirs(savepath)
    
        checkf = os.path.join(savepath,'cutout_'+str(sfid)+'.hdf5')
        npyf = os.path.join(savepath,'sub_'+str(sfid)+'.npy')
    else:
        checkf = 'cutout_'+str(sfid)+'.hdf5'
        npyf = 'sub_'+str(sfid)+'.npy'
    

    if os.path.lexists(checkf) and os.path.lexists(npyf) and clobber is False:
        if verbose:
            print("Subhalo cutout exists, skipping: ", checkf)
        download = False
        subobj = np.load(npyf)
        sub = subobj.all()   #?

        return checkf, sub, download
    else:
        if verbose:
            print("Getting subhalo cutout: ", checkf)

        #subhalo metadata
        try:
            sub = get(url)
            sim_obj=get(sim_url)
            snap_obj=get(snap_url)

            np.save(npyf,sub)
            file = get(url+"/cutout.hdf5",params,savepath)

            download = True
            #add attributes to header that Sunrise needs for GFM setting.. time, Omega, npart?
            with h5py.File(file,'a') as fo:
                header=fo['Header']
                header.attrs['Time']=1.0/(1.0 + snap_obj['redshift'])#actually scale factor
                header.attrs['HubbleParam']=sim_obj['hubble']
                header.attrs['Omega0']=sim_obj['omega_0']
                header.attrs['OmegaLambda']=sim_obj['omega_L']
                npart = [sub['len_gas'],sub['len_dm'],0,0,sub['len_stars'],sub['len_bhs']]
                header.attrs['NumPart_ThisFile']=np.asarray(npart)
                mtable=[0,sim_obj['mass_dm'],0,0,0,0]
                header.attrs['MassTable']=np.asarray(mtable)
                header.attrs['Redshift']=snap_obj['redshift']
        except:
            file = None
            sub = None
            download = False

        return file, sub, download


def process_subhalo(sim,snap,sfid,i,label,camera_obj,params=defaultparams,savepath=None,verbose=True,clobber=False,getlabel='StellarMass',resample=False):
    #defaultparams={'stars':'Coordinates,Velocities,GFM_StellarFormationTime,GFM_Metallicity,GFM_InitialMass,Masses'}

    file, sub, downloaded = get_subhalo(sim,snap,sfid,params=defaultparams,savepath=savepath,verbose=verbose,clobber=clobber)
    substuff = None

    if i % 100 ==0:
        print("Finished: ", i, snap, sfid, downloaded, (time.time() - start_time)/60.0, file)
        sys.stdout.flush()



    if file is None:
        print("Subhalo Failed: ", i, snap, sfid, downloaded, file)
        sys.stdout.flush()
        return (file, substuff)


    #check if resampled cutout exists, if so don't bother loading
    dirname = os.path.dirname(file)
    resampfile = os.path.join(dirname,'resampled_'+str(sfid)+'.hdf5')
    resamp_exists = os.path.lexists(resampfile)

    #if not resamp_exists:
        #load HDF5 file relevant fields
    #    with h5py.File(file) as h5f:
    #        masses = np.asarray(h5f['PartType4']['Masses'])*(1.0e10)/ilh  #convert to solar masses
    #        formfactors = np.asarray(h5f['PartType4']['GFM_StellarFormationTime'])  #in scale factor
    #        xpos = np.asarray(h5f['PartType4']['Coordinates'][:,0])-pos[0] #in ckpc/h
    #        ypos = np.asarray(h5f['PartType4']['Coordinates'][:,1])-pos[1]
    #        zpos = np.asarray(h5f['PartType4']['Coordinates'][:,2])-pos[2]

        #print masses.shape, ypos.shape, np.mean(ypos)


    #resample star particles
    

    #project into an image to return... or output Sunrise input files as needed
    #return sunrise directory... relative to savepath?
    sundir = os.path.join(dirname,label)
    if not os.path.lexists(sundir):
        os.mkdir(sundir)
        
    #sfrhist -- translate to center of subhalo

    #mcrx -- move cameras to negative subhalo position, set FOV smartly
    #camera position = -camdir*distance; distance from lightcone catalog Z and tz ? or camdir dot XYZ ?

    #broadband -- always redshift?
    

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

def get_lightcone_images_threaded(lcfile,geofile,sim='Illustris-2',clobber=False,savepath='/astro/snyder_lab2/Illustris',Np=2,maxq=10000,lim=None,label='FIELD'):
    data = ascii.read(lcfile)
    snapnums = np.asarray(data['col1'],dtype='str')
    sfids = np.asarray(data['col2'],dtype='str')
    x_cmpc = np.asarray(data['col13'])
    y_cmpc = np.asarray(data['col14'])
    z_cmpc = np.asarray(data['col15'])
    coords_ra = np.asarray(data['col3']) #degrees
    coords_dec = np.asarray(data['col4']) #degrees

    catalog = ilc.process_lightcone_catalog(lightcone=geofile,basedir='.')
    print(catalog.delb_arcmin)


    NUMBER_OF_PROCESSES=Np
    task_queue = Queue()
    done_queue = Queue()
    TASKS = []
    TASKS_DONE = []
    TASKS_LEFT = []

    N_objects = snapnums.shape[0]
    #figure out how to drain the queue before it fills the pipe

    if lim is None:
        lim=np.int64(N_objects)

    print("Items to process: ", lim)

    for i,sn in enumerate(snapnums[0:lim]):
        this_sfid = sfids[i]
        this_x = x_cmpc[i]
        this_y = y_cmpc[i]
        this_z = z_cmpc[i]
        this_ra = coords_ra[i]
        this_dec = coords_dec[i]

        task = (process_subhalo,(sim,sn,this_sfid,i,label,None))
        if i <= maxq:
            task_queue.put(task)
            TASKS.append(task)
        else:
            TASKS_LEFT.append(task)

    for p in range(NUMBER_OF_PROCESSES):
        Process(target=worker,args=(task_queue,done_queue),kwargs={'savepath':savepath,'verbose':False,'clobber':clobber}).start()

    cutout_files = []

    while len(TASKS_LEFT) > 0:
        cutout_files.append(done_queue.get())
        newtask = TASKS_LEFT.pop()
        task_queue.put(newtask)

    for i in range(min(maxq,lim)):
        cutout_files.append(done_queue.get())


    #build up images/analysis here!
    print(cutout_files[0:5])
    print(cutout_files[-5:])
    print(len(cutout_files))

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


def do_lightcone_images(savepath=None):

    if savepath is None:
        print("Requires savepath specification!")
        exit()

    lcfile = os.path.expandvars('$HOME/oasis_project/Lightcones/Illustris-2_RADEC_hudfwide_75Mpc_7_6_xyz_corners.txt')
    geofile = os.path.expandvars('$HOME/oasis_project/Lightcones/hudfwide_75Mpc_7_6_fixedh_xyz_NEW.txt')
    label = 'FIELDA'

    #note, for Gordon, stage intermediate lightcone data on local scratch?  Do this outside python:  FASTER

    #if not os.path.lexists(final_savepath):
    #    os.mkdir(final_savepath)

    #print "Staging files... ", final_savepath, intermediate_savepath
    #staget = time.time()
    #shutil.move(final_savepath,intermediate_savepath)
    #staget = time.time() - staget
    #print "Staging took: ", staget, ' seconds'

    #st1 = time.time()
    #get_lightcone_images(lcfile,geofile,sim='Illustris-2',clobber=False,savepath=savepath)
    #et1 = time.time()

    st = time.time()
    result = get_lightcone_images_threaded(lcfile,geofile,sim='Illustris-2',clobber=False,savepath=savepath,Np=8,lim=None,maxq=10000,label=label)
    et = time.time()

    #return files to temp project space
    #print "UnStaging files... ", intermediate_savepath, final_savepath
    #staget = time.time()
    #shutil.move(intermediate_savepath,final_savepath)
    #staget = time.time() - staget
    #print "Un-Staging took: ", staget, ' seconds'
    

    print('Threaded calculation took: ', et-st, ' seconds')
    #print 'Serial calculation took: ', et1-st1, ' seconds'


if __name__=="__main__":
    if len(sys.argv) != 2:
        print("Usage:  python illustris_api_utils_gordon.py SAVEPATH")
        exit()
    else:
        savepath=sys.argv[1]

    print("Saving intermediate files at: ", savepath)

    do_lightcone_images(savepath=savepath)
