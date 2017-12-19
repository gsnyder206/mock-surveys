import os
import sys
import numpy as np
import glob
import gfs_sublink_utils as gsu
import shutil
import math
import astropy
import astropy.io.fits as fits
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


sq_arcsec_per_sr = 42545170296.0
c = 3.0e8

lcfile_cols={'col1':'snapshot',
             'col2':'SubfindID',
             'col3':'ra_deg',
             'col4':'dec_deg',
             'col5':'ra_kpc',
             'col6':'dec_kpc',
             'col7':'ra_kpc_inferred',
             'col8':'dec_kpc_inferred',
             'col9':'true_z',
             'col10':'inferred_z',
             
             'col11':'peculiar_z',
             'col12':'true_kpc_per_arcsec',
             'col13':'X_cmpc',
             'col14':'Y_cmpc',
             'col15':'Z_cmpc',
             'col16':'ADD_cmpc',
             'col17':'ADD_cmpc_inferred',
             'col18':'snapshot_z',
             'col19':'geometric_z',
             'col20':'cylinder_number',

             'col21':'mstar_msun_rad',
             'col22':'mgas_msun_rad',
             'col23':'subhalo_mass_msun',
             'col24':'bhmass_msun_rad',
             'col25':'mbary_msun_rad',
             'col26':'sfr_msunperyr_rad',
             'col27':'bhrate_code',
             'col28':'camX_mpc',
             'col29':'camY_mpc',
             'col30':'camZ_mpc',
             
             'col31':'g_AB_absmag',
             'col32':'r_AB_absmag',
             'col33':'i_AB_absmag',
             'col34':'z_AB_absmag',
             'col35':'v_kms_camX',
             'col36':'v_kms_camY',
             'col37':'v_kms_camZ',
             'col38':'v_kms_hubble',
             'col39':'g_AB_appmag'}

scale_window=1.2
mass_window=2.0

def process_single_filter_subimage(image_parameters,galaxy_data,lcdata,filname,lambda_eff_microns):

    print('**** Processing subimage:  ',image_parameters)

    single_filter_subimage=np.ndarray((image_parameters['Npix'],image_parameters['Npix']))
    print(single_filter_subimage.shape)

    #find sources in this sub-image

    #use buffer to clear edge effects from lower left
    buf_deg=10.0/3600.0 #10 arcsec buffer?
    sub_indices=(lcdata['ra_deg']>=image_parameters['x1_deg']-buf_deg)*(lcdata['ra_deg']<image_parameters['x2_deg'])*(lcdata['dec_deg']>=image_parameters['y1_deg']-buf_deg)*(lcdata['dec_deg']<image_parameters['y2_deg'])

    sub_data=lcdata[sub_indices]

    success=0

    image_catalog={'filter':filname}

    #final source info
    xcen_list=[]
    ycen_list=[]
    final_flux_njy_list=[]
    ab_appmag_list=[]
    
    #galaxy_data=galaxy_data.fromkeys(['image_dir','scale','simlabel','Mvir','Mstar','Rhalf_stars'])

    #original image source data
    final_file_list=[]
    image_dir_list=[]
    scalefactor_list=[]
    simlabel_list=[]
    Mvir_list=[]
    mstar_list=[]
    rhalf_list=[]

    #lightcone entry data... all of it???
    
    
    for i,entry in enumerate(sub_data):
        #need pos_i, pos_j
        pos_i=np.int64( (entry['ra_deg']-image_parameters['x1_deg'])*np.float64(image_parameters['Npix'])/(image_parameters['x2_deg']-image_parameters['x1_deg'])   )
        pos_j=np.int64( (entry['dec_deg']-image_parameters['y1_deg'])*np.float64(image_parameters['Npix'])/(image_parameters['y2_deg']-image_parameters['y1_deg'])   )

        #select image file to insert
        mass_value=entry['subhalo_mass_msun']
        scale_value=1.0/(1.0+entry['true_z'])
        mstar_value=entry['mstar_msun_rad']
        

        galaxy_scale_indices=(galaxy_data['scale']>=scale_value/scale_window)
        galaxy_scale_indices*=(galaxy_data['scale']<=scale_value*scale_window)
        
        galaxy_mass_indices=(galaxy_data['Mvir']>=mass_value/mass_window)
        galaxy_mass_indices*=(galaxy_data['Mvir']<=mass_value*mass_window)

        galaxy_search_indices=galaxy_scale_indices*galaxy_mass_indices

        found_galaxy=False
        
        if np.sum(galaxy_search_indices)==0 and np.sum(galaxy_scale_indices)>0:
            #pick random galaxy and resize?
            #index into list of Trues
            random_index=random.randint(np.sum(galaxy_scale_indices))
            scale_where=np.where(galaxy_scale_indices==True)[0]
            galaxy_index=scale_where[random_index]
            success+=1
            found_galaxy=True
        elif np.sum(galaxy_search_indices)==0 and np.sum(galaxy_scale_indices)==0:
            galaxy_index=None
            pass
        else:
            random_index=random.randint(np.sum(galaxy_search_indices))
            galaxy_where=np.where(galaxy_search_indices==True)[0]
            galaxy_index=galaxy_where[random_index]
            success+=1
            found_galaxy=True


        
            
        found_data=False
        #now we have galaxy_index
        if galaxy_index is not None:
            mstar_factor=mstar_value/(galaxy_data['Mstar'][galaxy_index])
            size_factor=(mstar_factor)**0.5

            folder=galaxy_data['image_dir'][galaxy_index]
            label=galaxy_data['simlabel'][galaxy_index]
            possible_files=np.sort(np.asarray(glob.glob(folder+'/hires_images_cam??/'+label+'cam??_'+filname+'*.fits')))
            

            if possible_files.shape[0]>0:
                #pick a random camera
                file_index=random.randint(possible_files.shape[0])
                this_file=possible_files[file_index]
                try:
                    #locate file and load
                    found_data=True
                    this_hdu=fits.open(this_file)[0]
                    this_image=this_hdu.data
                    
                    pixstr=this_file.split('_')[-2][3:]

                    #this was foolish!

                    pixsize_arcsec=np.float64(pixstr)


                except:
                    found_data=False


                if found_data==True:
                    #these are in nJy-- preserve integral!

                    original_flux=np.sum(this_image)
                    
                    total_flux=original_flux*mstar_factor  #shrink luminosity by mstar factor

                    this_npix=this_image.shape[0]
                    
                    #resize--preserve proper units

                    desired_npix=np.int32( this_npix*(pixsize_arcsec/image_parameters['pix_size_arcsec'])*size_factor  )

                    resized_image=congrid.congrid(this_image,(desired_npix,desired_npix))
                    resized_flux=np.sum(resized_image)

                    resized_image=resized_image*(total_flux/resized_flux)


                    #is there a way to detect and avoid edge effects here??

                    #1. smooth image
                    #2. locate peak flux
                    #3. apply gaussian/exponential factor strong enough to eliminate to full-res image?
                    #4. use size info... to make sensible??

                    
                    #add to image

                    npsub=desired_npix
                    i1=pos_i
                    i2=pos_i+npsub
                    j1=pos_j
                    j2=pos_j+npsub
                    icen=np.float64(pos_i)+np.float64(npsub)/2.0
                    jcen=np.float64(pos_j)+np.float64(npsub)/2.0

                    #determine overlap image_parameters['Npix']
                    im0=0
                    im1=image_parameters['Npix']

                    sub_image_to_add=resized_image[im0-i1:npsub-(i2-im1),im0-j1:npsub-(j2-im1)]
                    print('Resized: ', mstar_factor, size_factor, pixsize_arcsec, this_npix, desired_npix, resized_image.shape,sub_image_to_add.shape, i1, i2, j1, j2)

                    if sub_image_to_add.shape[0] > 0 and sub_image_to_add.shape[1] > 0:
                        iplace=np.max([i1,0])
                        jplace=np.max([j1,0])
                        new_npi=sub_image_to_add.shape[0]
                        new_npj=sub_image_to_add.shape[1]
                        single_filter_subimage[iplace:iplace+new_npi,jplace:jplace+new_npj] += sub_image_to_add


                    if icen >= 0.0 and jcen >= 0.0:
                        #assemble catalog entries for this object
                        #final source info
                        xcen_list.append(icen)
                        ycen_list.append(jcen)
                        
                        final_flux_njy_list.append(total_flux)
                        ab_appmag_list.append(-2.5*np.log10((1.0e9)*(total_flux)/3632.0))
                        
                        #galaxy_data=galaxy_data.fromkeys(['image_dir','scale','simlabel','Mvir','Mstar','Rhalf_stars'])
                        
                        #original image source data
                        final_file_list=[]
                        image_dir_list=[]
                        scalefactor_list=[]
                        simlabel_list=[]
                        Mvir_list=[]
                        mstar_list=[]
                        rhalf_list=[]
                        
                        #lightcone entry data... all of it???
                            

                        
        print('found galaxy? ', found_galaxy, 'found data? ', found_data)
                


            
            
        
    print('**** Subimage successes: ', str(success), '  out of  ', i+1)

    
    return single_filter_subimage


def build_lightcone_images(lightcone_file,run_type='images',lim=None,minz=None,
                           image_filelabel='hlsp_misty_mockluvoir',
                           total_images=16,
                           do_sub_image_x=0,
                           do_sub_image_y=0,
                           n_cameras=19):

    #image parameters
    image_parameters={}
    image_parameters['Npix']=8192
    image_parameters['pix_size_arcsec']=0.00732421875
    #gives exactly 1-arcmin images

    lightcone_fov_arcmin=5.305160

    #give some buffer for clearing edge effects
    full_image_min=-0.95*lightcone_fov_arcmin/2.0
    full_image_max=0.95*lightcone_fov_arcmin/2.0

    x_min=full_image_min + (do_sub_image_x+0)*image_parameters['Npix']*image_parameters['pix_size_arcsec']/60.0
    x_max=full_image_min + (do_sub_image_x+1)*image_parameters['Npix']*image_parameters['pix_size_arcsec']/60.0

    y_min=full_image_min + (do_sub_image_y+0)*image_parameters['Npix']*image_parameters['pix_size_arcsec']/60.0
    y_max=full_image_min + (do_sub_image_y+1)*image_parameters['Npix']*image_parameters['pix_size_arcsec']/60.0


    image_parameters['x1_deg']=x_min/60.0
    image_parameters['x2_deg']=x_max/60.0
    image_parameters['y1_deg']=y_min/60.0
    image_parameters['y2_deg']=y_max/60.0
    image_parameters['xsub']=do_sub_image_x
    image_parameters['ysub']=do_sub_image_y

    #pick random seed 
    random.seed(image_parameters['xsub']+int(total_images**0.5)*image_parameters['ysub'])
    
    
    lightcone_dir=os.path.abspath(os.path.dirname(lightcone_file))
    print('Constructing lightcone data from: ', lightcone_file)

    output_dir = os.path.join(lightcone_dir,'luvoir_mosaics')
    print('Saving lightcone outputs in: ', output_dir)
    if not os.path.lexists(output_dir):
        os.mkdir(output_dir)
    
    lcdata=ascii.read(lightcone_file)

    for colkey in lcfile_cols:
        col_obj=lcdata[colkey]
        col_obj.name=lcfile_cols[colkey]
    
    filters = ['Z']#['U','B','V','I','Z','nircam_f115w','nircam_f150w','nircam_f200w']
    output_filters=['f850lp']#['f336w','f435w','f606w','f775w','f850lp','f115w','f150w','f200w']
    lambda_eff_microns = [0.35,0.45,0.55,0.75,0.85,1.15,1.50,2.0]


    #get galprops catalog(s)
    galaxy_data={}
    galaxy_data=galaxy_data.fromkeys(['image_dir','scale','simlabel','Mvir','Mstar','Rhalf_stars'])
    print(galaxy_data.keys())
    
    galprops_files=np.sort(np.asarray(glob.glob(lightcone_dir+'/galprops_data/VELA??_galprops.npy')))
    for gp_file in galprops_files:
        print('   ', gp_file)
        
        this_gp_dict=(np.load(gp_file,encoding='bytes')).all()

        snap_files=this_gp_dict[b'snap_files']
        sim_dir_list = np.str_(snap_files).split('/')[7::8]   #OMG wut
        sim_labels=[]
        sim_names=[]

        for sn in sim_dir_list:
            sim_labels.append(sn.rstrip('_sunrise'))
            sim_names.append(lightcone_dir+'/'+sn[0:6])
            
        sim_labels=np.asarray(sim_labels)
        sim_names=np.asarray(sim_names)
                
        assert(sim_labels.shape[0]==len(sim_names))

        if galaxy_data['image_dir'] is None:
            galaxy_data['image_dir']=sim_names
            galaxy_data['scale']=this_gp_dict[b'scale']
            galaxy_data['simlabel']=sim_labels
            galaxy_data['Mvir']=this_gp_dict[b'Mvir_dm']
            galaxy_data['Mstar']=this_gp_dict[b'stars_total_mass']
            galaxy_data['Rhalf_stars']=this_gp_dict[b'stars_rhalf']

        else:
            galaxy_data['image_dir']=np.append(galaxy_data['image_dir'],sim_names)
            galaxy_data['scale']=np.append(galaxy_data['scale'],this_gp_dict[b'scale'])
            galaxy_data['simlabel']=np.append(galaxy_data['simlabel'],sim_labels)
            galaxy_data['Mvir']=np.append(galaxy_data['Mvir'],this_gp_dict[b'Mvir_dm'])
            galaxy_data['Mstar']=np.append(galaxy_data['Mstar'],this_gp_dict[b'stars_total_mass'])
            galaxy_data['Rhalf_stars']=np.append(galaxy_data['Rhalf_stars'],this_gp_dict[b'stars_rhalf'])
                    
    assert(galaxy_data['image_dir'].shape==galaxy_data['Mvir'].shape)
    
               
    for i,filname in enumerate(filters):
        print('processing.. ', filname)
        
        single_filter_subimage=process_single_filter_subimage(image_parameters,galaxy_data,lcdata,filname,lambda_eff_microns[i])

        output_filename=os.path.join(output_dir,image_filelabel+'_imager'+'_fielda-subimage'+str(do_sub_image_x)+str(do_sub_image_y)+'-9-8_'+output_filters[i]+'_v1_lightcone.fits' )
        print('**** Saving Subimage: ', output_filename )

        outhdu=fits.PrimaryHDU(single_filter_subimage)
        outhdu.header['PIXSIZE']=(image_parameters['pix_size_arcsec'],'arcsec')
        
        outlist=fits.HDUList([outhdu])
        outlist.writeto(output_filename,overwrite=True)
    
        
    return
