import os
import sys
import numpy as np
import glob
import gfs_sublink_utils as gsu
import shutil
import math
import astropy
import astropy.io.fits as pyfits

#this function takes as inputs the return values of illustris_api_utils.get_subhalo()
#this is a snapshot file and subhalo metadata object


def setup_sunrise_lightcone(snap_cutout,subhalo_object,label,this_z,geofile,pos_mpc,submitcount,savepath,append=True,verbose=True,clobber=True,
                                        data_dir='$HOME/sunrise_data/',stub_dir='$HOME/PythonCode/mock-surveys/stubs_lightcones/',
                                        nthreads=24,use_scratch=True,run_type='images',pixsize_arcsec=0.032,rad_fact=10.0):

    #TBD--smart FOV and Npix customization-- want large enough field at ~JWST resolution?  0.016 arcsec?  half NC-short pix?
    # fov = 2x max extent?

    fits_file = os.path.abspath(snap_cutout)
    galprops_data = subhalo_object
    
    snap_dir = os.path.dirname(fits_file)
    
    stub_dir = os.path.expandvars(stub_dir)
    data_dir = os.path.expandvars(data_dir)
    
    print("Using stubs in.. ",stub_dir)
    stub_files = np.asarray(glob.glob(os.path.join(stub_dir,'*')))

    real_redshift=gsu.redshift_from_snapshot( subhalo_object['snap'] )
    scale_convert=(1.0/(gsu.ilh*(1.0 + real_redshift)))

    #put object at "observed" redshift
    redshift=this_z
    

    halfmassradstars = subhalo_object['halfmassrad_stars']/(gsu.ilh*(1.0 + real_redshift))   #physical kpc
    mstar = subhalo_object['mass_stars']*(1.0e10)/gsu.ilh

    used_rad_fact=min(100.0,rad_fact+mstar/1.0e9)  #scale with mass up to 100x halfmassrad

    desired_fov_kpc = min( used_rad_fact*halfmassradstars, 500.0 )  #max out at 0.5 Mpc

    desired_fov_arcsec = desired_fov_kpc/(gsu.illcos.kpc_proper_per_arcmin(this_z).value/60.0)
    npix_int = np.int64(desired_fov_arcsec/pixsize_arcsec)+1
    final_fov_arcsec = npix_int*pixsize_arcsec
    final_fov_kpc = final_fov_arcsec*(gsu.illcos.kpc_proper_per_arcmin(this_z).value/60.0)

    nrays=max( np.int64( (1.0e7)*((npix_int/200.0)**2.0)**0.67 ), np.int64(1e7) )


    mgas = subhalo_object['mass_gas']
    if mgas==0.0:
        aux_part_only=True
    else:
        aux_part_only=False


    nthreads=str(nthreads)

    run_dir = snap_dir+'/'+run_type+'_'+label

    sfid=subhalo_object['id']
        
    print('\tGenerating sunrise.sbatch file for %s...'%run_type)
    sbatch_fn   = 'sunrise_'+run_type+'_'+label+'_'+str(submitcount)+'.sbatch'
    
    final_fn = generate_sbatch_lightcone(run_dir = run_dir, snap_dir = snap_dir, filename = sbatch_fn, 
                                         galprops_data = galprops_data, run_type = run_type,savepath=savepath,ncpus=nthreads,walltime='08:00:00',use_scratch=use_scratch,append=append,isnap=sfid)



    if not os.path.lexists(run_dir):
        os.mkdir(run_dir)

    for sf in stub_files:
        shutil.copy(sf,run_dir)
            
    sfrhist_fn   = 'sfrhist.config'
    sfrhist_stub = 'sfrhist_base.stub'

    generate_sfrhist_config(run_dir = run_dir, filename = sfrhist_fn, data_dir=data_dir,
                            stub_name = sfrhist_stub,  fits_file = fits_file, 
                            galprops_data = galprops_data, run_type = run_type,
                            nthreads=nthreads,scale_convert=scale_convert,use_scratch=use_scratch,isnap=sfid)



    generate_campos(run_dir,this_z,geofile,pos_mpc,fov_kpc=final_fov_kpc)
    
    mcrx_fn   = 'mcrx.config'
    mcrx_stub = 'mcrx_base.stub'

    generate_mcrx_config(run_dir = run_dir, snap_dir = snap_dir, filename = mcrx_fn, 
                         stub_name = mcrx_stub,
                         galprops_data = galprops_data, run_type = run_type, nthreads=nthreads, 
                         cam_file='cam_pos.txt' ,use_scratch=use_scratch,isnap=sfid,npix=npix_int,aux_part_only=aux_part_only,nrays=nrays)


        
    if run_type == 'images': 
        broadband_fn   = 'broadband.config'
        broadband_stub = 'broadband_base.stub'

        generate_broadband_config_images(run_dir = run_dir, snap_dir = snap_dir, data_dir=data_dir, filename = broadband_fn, 
                                         stub_name = broadband_stub, 
                                         galprops_data = galprops_data,redshift=redshift,use_scratch=use_scratch,isnap=sfid)
    if run_type == 'grism': 
        broadband_fn   = 'broadband.config'
        broadband_stub = 'broadband_base.stub'

        generate_broadband_config_grism(run_dir = run_dir, snap_dir = snap_dir, data_dir=data_dir, filename = broadband_fn, 
                                        stub_name = broadband_stub, 
                                        galprops_data = galprops_data,redshift=redshift,use_scratch=use_scratch,isnap=sfid)






    
    #setup sfrhist
    #setup mcrx
    #setup broadband

    #TEST WITH SUBSET OF OBJECTS/CATALOGS
    
    return_dict={}
    return_dict['submitfile']=final_fn
    return_dict['this_npix']=npix_int
    return_dict['fov_kpc']=final_fov_kpc
    return_dict['run_dir']=run_dir
    return_dict['rad_fact']=used_rad_fact
    return_dict['nrays']=nrays

    return return_dict




def generate_campos(run_dir,this_z,geofile,pos_mpc,fov_kpc=120.0):

    lines = open(geofile,'r')
    for l in lines:
        if "Comoving Single Box L" in l:
            L_comoving = np.float32(l.split()[-1])
            L_comovingh = round(L_comoving*gsu.ilh,4)

        if "Delta Unit Vector" in l:
            ss = l.split("[")[-1].split("]")[0].split()
            xs = ss[0]
            ys = ss[1]
            zs = ss[2]   

        if "Direction Unit Vector" in l:
            ss = l.split("[")[-1].split("]")[0].split()
            xd = ss[0]
            yd = ss[1]
            zd = ss[2]
        if "del B" in l:
            delb_arcmin = np.float32(l.split()[-1])
        if "del A" in l:
            dela_arcmin = np.float32(l.split()[-1])
    lines.close()
    assert xs is not None
    assert xd is not None
    assert L_comoving is not None

    #camdir
    #just the direction unit vector from the lightcone file
    camdir_x = np.float32(xd)
    camdir_y = np.float32(yd)
    camdir_z = np.float32(zd)
    #camup
    #just the delta unit vector from the lightcone file
    camup_x = np.float32(xs)
    camup_y = np.float32(ys)
    camup_z = np.float32(zs)

    #I'm not positive what's the best option here?
    camdist_mpc = pos_mpc['z'] #(pos_mpc['x']**2+pos_mpc['y']**2+pos_mpc['z']**2)**0.5   #physical mpc in camera coords
    
    campos_x= -1.0*camdir_x*camdist_mpc*10.0**3  #want in kpc
    campos_y= -1.0*camdir_y*camdist_mpc*10.0**3
    campos_z= -1.0*camdir_z*camdist_mpc*10.0**3

    #update camdir here to dither image?
    
    fov_radians = math.asin(fov_kpc/(camdist_mpc*1.0e3))
    
    cf = open(run_dir+'/cam_pos.txt','w+')  #pos, dir, up, fov radians
    cf.write('{:12.4f} {:12.4f} {:12.4f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.10f}\n'.format(campos_x,campos_y,campos_z,camdir_x,camdir_y,camdir_z,camup_x,camup_y,camup_z,fov_radians))
    cf.write('\n')
    cf.close()
    
    return


def generate_sbatch_lightcone(run_dir, snap_dir, filename, galprops_data, run_type, savepath, ncpus='24', queue='compute',
                              email='gsnyder@stsci.edu',walltime='04:00:00',isnap=0,account='hsc102',use_scratch=False,append=True):

    
    filepath=savepath+'/'+filename

    if append is True:
        bsubf = open(filepath,'a')
    else:
        bsubf = open(filepath, 'w+')
        
        bsubf.write('#!/bin/bash\n')
        bsubf.write('#SBATCH -A '+account+'\n')
        bsubf.write('#SBATCH --partition='+queue+'\n')
        bsubf.write('#SBATCH -t '+walltime+'\n')
        bsubf.write('#SBATCH --nodes=1\n')
        bsubf.write('#SBATCH --ntasks-per-node='+ncpus+'\n')
    
        bsubf.write('#SBATCH --export=ALL\n')
        bsubf.write('#SBATCH --job-name=sunrise_'+run_type+'\n')
        bsubf.write('#SBATCH --output='+filepath.rstrip('.sbatch')+'_slurm.out'+'\n')
        bsubf.write('\n')
    
    bsubf.write('cd '+run_dir+' \n')   #go to directory where job should run
    if use_scratch is True:
        bsubf.write('mkdir -p /scratch/$USER/$SLURM_JOBID/'+str(isnap)+'\n')
    bsubf.write('/home/gsnyder/bin/sfrhist sfrhist.config > sfrhist.out 2> sfrhist.err\n')
    bsubf.write('/home/gsnyder/bin/mcrx mcrx.config > mcrx.out 2> mcrx.err\n')
    if run_type=='images':
        bsubf.write('/home/gsnyder/bin/broadband broadbandz.config > broadbandz.out 2> broadbandz.err\n')
        bsubf.write('/home/gsnyder/bin/broadband broadband.config > broadband.out 2> broadband.err\n')
        if use_scratch is True:
            bsubf.write('rm -rf /scratch/$USER/$SLURM_JOBID/'+str(isnap)+'/sfrhist.fits\n')   #enable this after testing
            bsubf.write('rm -rf /scratch/$USER/$SLURM_JOBID/'+str(isnap)+'/mcrx.fits\n')   #enable this after testing
        #bsubf.write(os.path.expandvars('python $SYNIMAGE_CODE/candelize.py\n'))
    elif run_type=='ifu':
        bsubf.write('rm -rf sfrhist.fits\n')   #enable this after testing
    elif run_type=='grism':
        bsubf.write('/home/gsnyder/bin/broadband broadbandgrism.config > broadbandgrism.out 2> broadbandgrism.err\n')

    if use_scratch is True:
        bsubf.write('cp /scratch/$USER/$SLURM_JOBID/'+str(isnap)+'/*.fits .\n')

    bsubf.write('tar cf textfiles.tar simpar *.config *.stub *.txt *.out *.err *.hpmout --remove-files\n')

    bsubf.write('\n')
    bsubf.write('\n')
    bsubf.close()

    return filepath





def setup_sunrise_illustris_subhalo(snap_cutout,subhalo_object,verbose=True,clobber=True,
                                    stub_dir='$HOME/PythonCode/mock-surveys/stubs_illustris/',
                                    data_dir='$HOME/sunrise_data/',
                                    nthreads=24,redshift_override=None,walltime_limit='02:00:00',use_scratch=True):

    fits_file = os.path.abspath(snap_cutout)
    galprops_data = subhalo_object
    
    snap_dir = os.path.dirname(fits_file)
    
    print("Setting up sunrise run in.. ", snap_dir)

    stub_dir = os.path.expandvars(stub_dir)
    data_dir = os.path.expandvars(data_dir)
    
    print("Using stubs in.. ",stub_dir)
    stub_files = np.asarray(glob.glob(os.path.join('stub_dir','*')))

    
    list_of_types = ['grism','images']

    idx=None

    real_redshift=gsu.redshift_from_snapshot( subhalo_object['snap'] )
    scale_convert=(1.0/(gsu.ilh*(1.0 + real_redshift)))
    
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

        generate_sfrhist_config(run_dir = run_dir, filename = sfrhist_fn, data_dir=data_dir,
                                stub_name = sfrhist_stub,  fits_file = fits_file, 
                                galprops_data = galprops_data, run_type = run_type,
                                nthreads=nthreads, idx = idx,scale_convert=scale_convert,use_scratch=use_scratch)


        print('\tGenerating mcrx.config file for %s...'%run_type)
        mcrx_fn   = 'mcrx.config'
        mcrx_stub = os.path.join(stub_dir,'mcrx_base.stub')

        generate_mcrx_config(run_dir = run_dir, snap_dir = snap_dir, filename = mcrx_fn, 
                             stub_name = mcrx_stub,
                             galprops_data = galprops_data, run_type = run_type, nthreads=nthreads, cam_file=None , idx = idx,use_scratch=use_scratch)


        
        if run_type == 'images': 
            print('\tGenerating broadband.config file for %s...'%run_type)
            broadband_fn   = 'broadband.config'
            broadband_stub = os.path.join(stub_dir,'broadband_base.stub')

            generate_broadband_config_images(run_dir = run_dir, snap_dir = snap_dir, data_dir=data_dir, filename = broadband_fn, 
                                             stub_name = broadband_stub, 
                                             galprops_data = galprops_data, idx = idx,redshift=redshift,use_scratch=use_scratch)
        if run_type == 'grism': 
            print('\tGenerating broadband.config file for %s...'%run_type)
            broadband_fn   = 'broadband.config'
            broadband_stub = os.path.join(stub_dir, 'broadband_base.stub')

            generate_broadband_config_grism(run_dir = run_dir, snap_dir = snap_dir, data_dir=data_dir, filename = broadband_fn, 
                                            stub_name = broadband_stub, 
                                            galprops_data = galprops_data, idx = idx,redshift=redshift,use_scratch=use_scratch)





        print('\tGenerating sunrise.sbatch file for %s...'%run_type)
        sbatch_fn   = 'sunrise.sbatch'		
        final_fn = generate_sbatch(run_dir = run_dir, snap_dir = snap_dir, filename = sbatch_fn, 
                                 galprops_data = galprops_data, run_type = run_type,ncpus=nthreads,walltime=walltime_limit,use_scratch=use_scratch)

        
    
    return





def setup_sunrise_enzo(snap_fits,prop_file,verbose=True,clobber=True,
                       stub_dir='$HOME/PythonCode/mock-surveys/stubs_enzo/',
                       data_dir='$HOME/sunrise_data/',
                       nthreads=24,redshift_override=None,walltime_limit='02:00:00',use_scratch=True):

    fits_file = os.path.abspath(snap_fits)
    
    snap_dir = os.path.dirname(fits_file)
    
    print("Setting up sunrise run in.. ", snap_dir)

    stub_dir = os.path.expandvars(stub_dir)
    data_dir = os.path.expandvars(data_dir)
    
    print("Using stubs in.. ",stub_dir)
    stub_files = np.asarray(glob.glob(os.path.join(stub_dir,'*')))



    
    list_of_types = ['images','grism']

    galprops_data = np.load(prop_file)[()]
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

        generate_sfrhist_config(run_dir = run_dir, filename = sfrhist_fn, data_dir=data_dir,
                                stub_name = sfrhist_stub,  fits_file = fits_file, 
                                galprops_data = galprops_data, run_type = run_type,
                                nthreads=nthreads, idx = idx,scale_convert=scale_convert,use_scratch=use_scratch,use_cf00=True)


        print('\tGenerating mcrx.config file for %s...'%run_type)
        mcrx_fn   = 'mcrx.config'
        mcrx_stub = os.path.join(stub_dir,'mcrx_base.stub')

        generate_mcrx_config(run_dir = run_dir, snap_dir = snap_dir, filename = mcrx_fn, 
                             stub_name = mcrx_stub,npix=1000,
                             galprops_data = galprops_data, run_type = run_type, nthreads=nthreads, cam_file=cam_file , idx = idx,use_scratch=use_scratch, aux_part_only=True)


        
        if run_type == 'images': 
            print('\tGenerating broadband.config file for %s...'%run_type)
            broadband_fn   = 'broadband.config'
            broadband_stub = os.path.join(stub_dir,'broadband_base.stub')

            generate_broadband_config_images(run_dir = run_dir, snap_dir = snap_dir, data_dir=data_dir, filename = broadband_fn, 
                                             stub_name = broadband_stub, 
                                             galprops_data = galprops_data, idx = idx,redshift=redshift,use_scratch=use_scratch)
        if run_type == 'grism': 
            print('\tGenerating broadband.config file for %s...'%run_type)
            broadband_fn   = 'broadband.config'
            broadband_stub = os.path.join(stub_dir, 'broadband_base.stub')

            generate_broadband_config_grism(run_dir = run_dir, snap_dir = snap_dir, data_dir=data_dir, filename = broadband_fn, 
                                            stub_name = broadband_stub, 
                                            galprops_data = galprops_data, idx = idx,redshift=redshift,use_scratch=use_scratch)





        print('\tGenerating sunrise.sbatch file for %s...'%run_type)
        sbatch_fn   = 'sunrise.sbatch'		
        final_fn = generate_sbatch(run_dir = run_dir, snap_dir = snap_dir, filename = sbatch_fn, 
                                   galprops_data = galprops_data, run_type = run_type,ncpus=nthreads,walltime=walltime_limit,use_scratch=use_scratch,candelize=True)

        
    
    return




def generate_sfrhist_config(run_dir, filename, data_dir, stub_name, fits_file, galprops_data, run_type, nthreads='1', idx = None,scale_convert=1.0,use_scratch=False,isnap=None,use_cf00=False):
    if use_scratch is True:
        if isnap is not None:
            int_dir='/scratch/$USER/$SLURM_JOBID/'+str(isnap)
        else:
            int_dir='/scratch/$USER/$SLURM_JOBID'
    else:
        int_dir=run_dir

    sf = open(run_dir+'/'+filename,'w+')
    sf.write('#Parameter File for Sunrise, sfrhist\n\n')
    sf.write('include_file        		%s\n\n'%stub_name)
    sf.write('snapshot_file       		%s\n'%fits_file)
    sf.write('output_file          		%s\n\n'%(int_dir+'/sfrhist.fits'))
    sf.write('n_threads          		'+nthreads+'\n')
    
    #only one of these should be translated
    gridw=200
    if idx==None:
        sf.write('translate_origin          %.2f\t%.2f\t%.2f         / [kpc]\n'%(galprops_data['pos_x']*scale_convert, galprops_data['pos_y']*scale_convert, galprops_data['pos_z']*scale_convert))
        sf.write('grid_min			%.1f\t%.1f\t%.1f         / [kpc]\n'%(-1.0*gridw, -1.0*gridw,-1.0*gridw))
        sf.write('grid_max			%.1f\t%.1f\t%.1f         / [kpc]\n\n\n'%(1.0*gridw,1.0*gridw,1.0*gridw))
    else:
        sf.write('translate_origin          %.2f\t%.2f\t%.2f         / [kpc]\n'%(galprops_data['stars_maxndens'][idx][0], galprops_data['stars_maxndens'][idx][1], galprops_data['stars_maxndens'][idx][2]))




    if run_type == 'images':
        sf.write('min_wavelength			%s\n'%("0.02e-6"))
        sf.write('max_wavelength			%s\n\n'%("5.0e-6"))
        
        sf.write('mappings_sed_file			%s\n'%(data_dir+"Smodel-lores128.fits"))
        if use_cf00==True:
            sf.write('stellarmodelfile                     %s\n'%(data_dir+"bc03_cf00_1e7_new_ism_both.fits"))
        else:
            sf.write('stellarmodelfile			%s\n'%(data_dir+"Patrik-imfKroupa-Zmulti-ml.fits"))


    elif run_type == 'ifu':
        sf.write('min_wavelength			%s\n'%("0.6450e-6"))
        sf.write('max_wavelength			%s\n\n'%("0.6650e-6"))
        
        sf.write('mappings_sed_file			%s\n'%(data_dir+"Smodel_full_hires.fits"))
        sf.write('stellarmodelfile			%s\n'%(data_dir+"logspace-Patrik-imfKroupa-geneva-Zmulti-hires.fits"))

    elif run_type == 'grism':
        sf.write('min_wavelength			%s\n'%("4.5e-7"))  #NIRISS F200W
        sf.write('max_wavelength			%s\n\n'%("5.5e-7"))
        
        #sf.write('mappings_sed_file			%s\n'%(data_dir+'Mappings_Smodels_gfs.fits'))
        #sf.write('stellarmodelfile			%s\n'%(data_dir+'GFS_combined_nolines.fits'))   
        sf.write('mappings_sed_file			%s\n'%(data_dir+"Smodel_full_hires.fits"))
        sf.write('stellarmodelfile			%s\n'%(data_dir+"logspace-Patrik-imfKroupa-geneva-Zmulti-hires.fits"))



    sf.close()
    print('\t\tSuccessfully generated %s'%filename)
    
    return




def generate_mcrx_config(run_dir, snap_dir, filename, stub_name, galprops_data, run_type, nthreads='1',cam_file=None, idx = None,use_scratch=False,isnap=None,npix=400,aux_part_only=False,nrays=1e7):
    if use_scratch is True:
        if isnap is not None:
            int_dir='/scratch/$USER/$SLURM_JOBID/'+str(isnap)
        else:
            int_dir='/scratch/$USER/$SLURM_JOBID'
    else:
        int_dir=run_dir


    mf = open(run_dir+'/'+filename,'w+')



    mf.write('#Parameter File for Sunrise, mcrx\n\n')
    mf.write('include_file         %s\n\n'%stub_name)
    mf.write('input_file           %s\n'%(int_dir+'/sfrhist.fits'))
    mf.write('output_file          %s\n'%(int_dir+'/mcrx.fits'))

    #approximating Torrey and HST13887 settings
    if cam_file is None:
        mf.write('exclude_south_pole       true #comment\n')
        mf.write('camerafov     120\n')
        mf.write('ntheta        2\n')
        mf.write('nphi        3\n')
    else:
        mf.write('camera_positions      %s\n'%(cam_file))


    mf.write('n_threads          		'+nthreads+'\n')
        
    if run_type != 'ifu':
        mf.write('use_kinematics	   %s\n'%('false #True for IFU'))
    else:
        mf.write('use_kinematics	   %s\n'%('true #False for images'))

    #move npixels to .config file

    if run_type == 'images':
        mf.write('npixels     '+str(npix)+'\n')
    elif run_type == 'ifu':
        mf.write('npixels     '+str(npix)+'\n')
    else:
        mf.write('npixels     '+str(npix)+'\n')

    if aux_part_only is True:
        mf.write('aux_particles_only      true\n')
    else:
        mf.write('aux_particles_only      false\n')

    mf.write('nrays_nonscatter        '+str(nrays)+'\n')

    mf.close()



    print('\t\tSuccessfully generated %s'%filename)

    return


def generate_broadband_config_images(run_dir, snap_dir, data_dir, filename, stub_name, galprops_data, idx = None,redshift=0.0,use_scratch=False,isnap=None):

    #copy sunrise filter folder to snap_dir+'/inputs/sunrise_filters/'

    bf = open(run_dir+'/'+filename,'w+')
    red0=0.0
        

    if use_scratch is True:
        if isnap is not None:
            int_dir='/scratch/$USER/$SLURM_JOBID/'+str(isnap)
        else:
            int_dir='/scratch/$USER/$SLURM_JOBID'
    else:
        int_dir=run_dir


    bf.write('#Parameter File for Sunrise, broadband\n\n')
    bf.write('include_file                      %s\n\n'%stub_name)
    bf.write('redshift                          %.3f\n\n'%red0)
    bf.write('input_file                        %s\n'%(int_dir+'/mcrx.fits'))
    bf.write('output_file                       %s\n'%(int_dir+'/broadband.fits'))
    bf.write('filter_list                       %s\n'%(data_dir+'sunrise_filters/filters_rest'))
    bf.write('filter_file_directory             %s\n'%(data_dir+'sunrise_filters/'))
    bf.close()
    
    bfz = open(run_dir+'/'+filename.replace('broadband','broadbandz'),'w+')
    
    bfz.write('#Parameter File for Sunrise, broadband\n\n')
    bfz.write('include_file                      %s\n\n'%stub_name)
    bfz.write('redshift                          %.8f\n\n'%redshift)
    bfz.write('input_file                        %s\n'%(int_dir+'/mcrx.fits'))
    bfz.write('output_file                       %s\n'%(int_dir+'/broadbandz.fits'))
    bfz.write('filter_list                       %s\n'%(data_dir+'sunrise_filters/filters_st'))
    bfz.write('filter_file_directory             %s\n'%(data_dir+'sunrise_filters/'))
    bfz.close()
    



    print('\t\tSuccessfully generated %s'%filename)

    return


def generate_broadband_config_grism(run_dir, snap_dir, data_dir, filename, stub_name, galprops_data, idx = None, redshift=0.0,use_scratch=False,isnap=None):


    bfg = open(run_dir+'/'+filename.replace('broadband','broadbandgrism'),'w+')
    
    if use_scratch is True:
        if isnap is not None:
            int_dir='/scratch/$USER/$SLURM_JOBID/'+str(isnap)
        else:
            int_dir='/scratch/$USER/$SLURM_JOBID'
    else:
        int_dir=run_dir

    bfg.write('#Parameter File for Sunrise, broadband\n\n')
    bfg.write('include_file                      %s\n\n'%stub_name)
    bfg.write('redshift                          %.3f\n\n'%redshift)
    bfg.write('input_file                        %s\n'%(int_dir+'/mcrx.fits'))
    bfg.write('output_file                       %s\n'%(int_dir+'/grism.fits'))
    bfg.write('filter_list                       %s\n'%(data_dir+'sunrise_filters/filters_grism'))
    bfg.write('filter_file_directory             %s\n'%(data_dir+'sunrise_filters/'))
    bfg.close()
    
    
    
    
    print('\t\tSuccessfully generated %s'%filename)
    
    return


def generate_sbatch(run_dir, snap_dir, filename, galprops_data, run_type, ncpus='24', queue='compute',email='gsnyder@stsci.edu',walltime='04:00:00',isnap=0,account='hsc102',use_scratch=False, candelize=False):

    bsubf = open(run_dir+'/'+filename, 'w+')
    bsubf.write('#!/bin/bash\n')
    bsubf.write('#SBATCH -A '+account+'\n')
    bsubf.write('#SBATCH --partition='+queue+'\n')
    bsubf.write('#SBATCH -t '+walltime+'\n')
    bsubf.write('#SBATCH --nodes=1\n')
    bsubf.write('#SBATCH --ntasks-per-node='+ncpus+'\n')
    
    bsubf.write('#SBATCH --export=ALL\n')
    bsubf.write('#SBATCH --job-name=sunrise_'+run_type+'\n')
    bsubf.write('#SBATCH --output='+run_dir+'/sunrise_slurm.out\n')
    bsubf.write('\n')
    
    bsubf.write('cd '+run_dir+' \n')   #go to directory where job should run
    bsubf.write('/home/gsnyder/bin/sfrhist sfrhist.config > sfrhist.out 2> sfrhist.err\n')
    bsubf.write('/home/gsnyder/bin/mcrx mcrx.config > mcrx.out 2> mcrx.err\n')
    if run_type=='images':
        #for these, may want to use:  https://github.com/gsnyder206/synthetic-image-morph/blob/master/tng/filters_lsst_light.txt
        bsubf.write('/home/gsnyder/bin/broadband broadbandz.config > broadbandz.out 2> broadbandz.err\n')
        bsubf.write('/home/gsnyder/bin/broadband broadband.config > broadband.out 2> broadband.err\n')
        
        #bsubf.write(os.path.expandvars('python $SYNIMAGE_CODE/mock_panstarrs.py\n'))  #doesn't exist yet!
        #follow this code below to write mock_panstarrs:
        if candelize==True:
            bsubf.write(os.path.expandvars('python $SYNIMAGE_CODE/candelize_enzo.py\n'))
    elif run_type=='ifu':
        bsubf.write('rm -rf sfrhist.fits\n')   #enable this after testing
        #bsubf.write('gzip -9 mcrx.fits\n')
    elif run_type=='grism':
        bsubf.write('/home/gsnyder/bin/broadband broadbandgrism.config > broadbandgrism.out 2> broadbandgrism.err\n')
        bsubf.write('rm -rf sfrhist.fits\n')   #enable this after testing
        #bsubf.write('rm -rf mcrx.fits\n')   #enable this after testing

    if use_scratch is True:
        bsubf.write('cp /scratch/$USER/$SLURM_JOBID/broadband*.fits .')

    
    bsubf.write('\n')
    bsubf.close()

    return os.path.abspath(run_dir+'/'+filename)







#will need a generate_sbatch for slurm jobs on Comet
def generate_qsub(run_dir, snap_dir, filename, galprops_data, run_type, ncpus='12', model='wes', queue='normal',email='rsimons@jhu.edu',walltime='04:00:00',isnap=0):

    bsubf = open(run_dir+'/'+filename, 'w+')
    bsubf.write('#!/bin/bash\n')
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
    bsubf.write('/u/gfsnyder/bin/sfrhist sfrhist.config > sfrhist.out 2> sfrhist.err\n')
    bsubf.write('/u/gfsnyder/bin/mcrx mcrx.config > mcrx.out 2> mcrx.err\n')
    if run_type=='images':
        bsubf.write('/u/gfsnyder/bin/broadband broadbandz.config > broadbandz.out 2> broadbandz.err\n')
        bsubf.write('/u/gfsnyder/bin/broadband broadband.config > broadband.out 2> broadband.err\n')
        bsubf.write('rm -rf sfrhist.fits\n')   #enable this after testing
        bsubf.write('rm -rf mcrx.fits\n')   #enable this after testing
        #bsubf.write(os.path.expandvars('python $SYNIMAGE_CODE/candelize.py\n'))
        bsubf.write('pigz -9 -p '+str(ncpus)+' broadband.fits\n')
    elif run_type=='ifu':
        #bsubf.write('rm -rf sfrhist.fits\n')   #enable this after testing
        bsubf.write('gzip -9 mcrx.fits\n')
    elif run_type=='grism':
        bsubf.write('/u/gfsnyder/bin/broadband broadbandgrism.config > broadbandgrism.out 2> broadbandgrism.err\n')
        #bsubf.write('rm -rf sfrhist.fits\n')   #enable this after testing
        #bsubf.write('rm -rf mcrx.fits\n')   #enable this after testing
        
    bsubf.close()




    return os.path.abspath(run_dir+'/'+filename)
