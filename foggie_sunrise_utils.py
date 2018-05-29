import os
import sys
import numpy as np
import glob
import shutil
import math
import astropy
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
import illustris_sunrise_utils as isu

def setup_sunrise_enzo(snap_fits,prop_file,verbose=True,clobber=True,
                       stub_dir='$HOME/PythonCode/mock-surveys/stubs_enzo/',
                       data_dir='$HOME/sunrise_data/',
                       pleiades=True,nthreads=20,redshift_override=None,walltime_limit='02:00:00',use_scratch=False):

    fits_file = os.path.abspath(snap_fits)
    
    snap_dir = os.path.dirname(fits_file)
    
    print("Setting up sunrise run in.. ", snap_dir)

    stub_dir = os.path.expandvars(stub_dir)
    data_dir = os.path.expandvars(data_dir)
    
    print("Using stubs in.. ",stub_dir)
    stub_files = np.asarray(glob.glob(os.path.join(stub_dir,'*')))



    
    list_of_types = ['images','grism']

    galprops_data = np.load(prop_file)[()]
    print(galprops_data)
    
    #oh wow
    dirdir=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(fits_file))))
    dirdirdir=os.path.join(dirdir,os.path.basename(dirdir))
    print('checking.. ', dirdir, dirdirdir, fits_file)

    idx = np.argwhere( galprops_data['snap_files']==dirdirdir )[0][0]

    cam_file = fits_file.rstrip('.fits')+'.cameras'

    #real_redshift=gsu.redshift_from_snapshot( subhalo_object['snap'] )
    #scale_convert=(1.0/(gsu.ilh*(1.0 + real_redshift)))

    #I think this always works in physical kpc anyway?  think further if need other redshifts
    scale_convert=1.0

    real_redshift=pyfits.open(snap_fits)['YT'].header['REDSHIFT']

    if redshift_override is None:
        redshift=real_redshift
    else:
        redshift=redshift_override


    print('redshift= ',redshift)
    
    nthreads=str(nthreads)
        
    for run_type in list_of_types:
        run_dir = snap_dir+'/%s'%run_type
        if not os.path.lexists(run_dir):
            os.mkdir(run_dir)

        for sf in stub_files:
            shutil.copy(sf,run_dir)
            
        print('\tGenerating sfrhist.config file for %s...'%run_type)
        sfrhist_fn   = 'sfrhist.config'
        sfrhist_stub = os.path.join(stub_dir,'sfrhist_base.stub')

        isu.generate_sfrhist_config(run_dir = run_dir, filename = sfrhist_fn, data_dir=data_dir,
                                stub_name = sfrhist_stub,  fits_file = fits_file, 
                                galprops_data = galprops_data, run_type = run_type,
                                nthreads=nthreads, idx = idx,scale_convert=scale_convert,use_scratch=use_scratch,use_cf00=True)


        print('\tGenerating mcrx.config file for %s...'%run_type)
        mcrx_fn   = 'mcrx.config'
        mcrx_stub = os.path.join(stub_dir,'mcrx_base.stub')

        isu.generate_mcrx_config(run_dir = run_dir, snap_dir = snap_dir, filename = mcrx_fn, 
                             stub_name = mcrx_stub,npix=1000,
                             galprops_data = galprops_data, run_type = run_type, nthreads=nthreads, cam_file=cam_file , idx = idx,use_scratch=use_scratch, aux_part_only=True)


        
        if run_type == 'images': 
            print('\tGenerating broadband.config file for %s...'%run_type)
            broadband_fn   = 'broadband.config'
            broadband_stub = os.path.join(stub_dir,'broadband_base.stub')

            isu.generate_broadband_config_images(run_dir = run_dir, snap_dir = snap_dir, data_dir=data_dir, filename = broadband_fn, 
                                             stub_name = broadband_stub, 
                                             galprops_data = galprops_data, idx = idx,redshift=redshift,use_scratch=use_scratch)
        if run_type == 'grism': 
            print('\tGenerating broadband.config file for %s...'%run_type)
            broadband_fn   = 'broadband.config'
            broadband_stub = os.path.join(stub_dir, 'broadband_base.stub')

            isu.generate_broadband_config_grism(run_dir = run_dir, snap_dir = snap_dir, data_dir=data_dir, filename = broadband_fn, 
                                            stub_name = broadband_stub, 
                                            galprops_data = galprops_data, idx = idx,redshift=redshift,use_scratch=use_scratch)





        print('\tGenerating sunrise.sbatch file for %s...'%run_type)
        if pleiades is True:
            qsub_fn = 'sunrise.qsub'
            final_fn = generate_qsub(run_dir = run_dir, snap_dir = snap_dir, filename = qsub_fn, model='ivy',
                                     galprops_data = galprops_data, run_type = run_type,ncpus=nthreads,walltime=walltime_limit)
        else:
            sbatch_fn   = 'sunrise.sbatch'		
            final_fn = isu.generate_sbatch(run_dir = run_dir, snap_dir = snap_dir, filename = sbatch_fn, 
                                       galprops_data = galprops_data, run_type = run_type,ncpus=nthreads,walltime=walltime_limit,use_scratch=use_scratch,candelize=True)

        
    
    return


def generate_qsub(run_dir, snap_dir, filename, galprops_data, run_type, group='s1938',ncpus='12', model='ivy', queue='normal',email='gsnyder@stsci.edu',walltime='04:00:00',isnap=0):

    bsubf = open(run_dir+'/'+filename, 'w+')
    bsubf.write('#!/bin/bash\n')
    bsubf.write('#PBS -W group_list='+group+'\n')
    bsubf.write('#PBS -S /bin/bash\n')   #apparently this is a thing
    bsubf.write('#PBS -l select=1:ncpus='+ncpus+':model='+model+'\n')   #selects cpu model and number (sunrise uses 1 node)
    bsubf.write('#PBS -l walltime='+walltime+'\n')    #hh:mm:ss before job is killed
    bsubf.write('#PBS -q '+queue+'\n')       #selects queue to submit to 
    bsubf.write('#PBS -N sunrise_'+run_type+'\n')     #selects job name
    bsubf.write('#PBS -M '+email+'\n')  #notifies job info to this email address 
    bsubf.write('#PBS -m abe\n')  #set notification types (abe=abort, begin, end)
    bsubf.write('#PBS -o '+run_dir+'/sunrise_pbs.out\n')  #save standard output here
    bsubf.write('#PBS -e '+run_dir+'/sunrise_pbs.err\n')  #save standard error here
    bsubf.write('#PBS -V\n')    #export environment variables at start of job
    
    bsubf.write('cd '+run_dir+' \n')   #go to directory where job should run
    bsubf.write(os.path.expandvars('$HOME')+'/bin/sfrhist sfrhist.config > sfrhist.out 2> sfrhist.err\n')
    bsubf.write(os.path.expandvars('$HOME')+'/bin/mcrx mcrx.config > mcrx.out 2> mcrx.err\n')
    if run_type=='images':
        bsubf.write(os.path.expandvars('$HOME')+'/bin/broadband broadbandz.config > broadbandz.out 2> broadbandz.err\n')
        bsubf.write(os.path.expandvars('$HOME')+'/bin/broadband broadband.config > broadband.out 2> broadband.err\n')
        #bsubf.write('rm -rf sfrhist.fits\n')   #enable this after testing
        #bsubf.write('rm -rf mcrx.fits\n')   #enable this after testing
        #bsubf.write(os.path.expandvars('python $SYNIMAGE_CODE/candelize.py\n'))
        #bsubf.write('pigz -9 -p '+str(ncpus)+' broadband.fits\n')
    elif run_type=='ifu':
        #bsubf.write('rm -rf sfrhist.fits\n')   #enable this after testing
        bsubf.write('gzip -9 mcrx.fits\n')
    elif run_type=='grism':
        bsubf.write(os.path.expandvars('$HOME')+'/bin/broadband broadbandgrism.config > broadbandgrism.out 2> broadbandgrism.err\n')
        #bsubf.write('rm -rf sfrhist.fits\n')   #enable this after testing
        #bsubf.write('rm -rf mcrx.fits\n')   #enable this after testing
        
    bsubf.close()




    return os.path.abspath(run_dir+'/'+filename)


def setup_foggie_image_pipeline():

    dat=ascii.read('sunrise_enzo_snap_list.txt')
    snaps=dat['snapshot_list']
    props=dat['galprops_list']
    
    for snap,prop in zip(snaps,props):
        print(snap, prop)
        
        setup_sunrise_enzo(snap,prop)
    
    #for loop over all snapshots in directory

    return



if __name__=="__main__":
    setup_foggie_image_pipeline()
