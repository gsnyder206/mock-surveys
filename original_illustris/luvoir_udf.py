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


def apply_adaptive_smoothing():
    kfile='hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f200w_v1_lightcone.fits'


    adaptive_smooth('hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f336w_v1_lightcone.fits',kfile)
    adaptive_smooth('hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f435w_v1_lightcone.fits',kfile)
    adaptive_smooth('hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f606w_v1_lightcone.fits',kfile)
    adaptive_smooth('hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f775w_v1_lightcone.fits',kfile)
    adaptive_smooth('hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f850lp_v1_lightcone.fits',kfile)

    adaptive_smooth('hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f115w_v1_lightcone.fits',kfile)
    adaptive_smooth('hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f150w_v1_lightcone.fits',kfile)
    adaptive_smooth('hlsp_misty_mockluvoir_imager_fielda-subimage00-9-8_f200w_v1_lightcone.fits',kfile)

    return

def adaptive_smooth(targetfile,distancefile):

    fitsfile=os.path.abspath(targetfile)
    distancefile=os.path.abspath(distancefile)
    
    filterstr=fitsfile.split('_')[-3]
    
    new_fn=fitsfile.rstrip('_v1_lightcone.fits').rstrip(filterstr).rstrip('_')+'-smooth'+'_'+filterstr+'_v1_lightcone.fits'

    print('smoothing..', new_fn)

    hdulist=fits.open(fitsfile)
    header=hdulist[0].header
    data=np.float32(hdulist[0].data)

    distlist=fits.open(distancefile)
    distdata=np.float32(distlist[0].data)

    #measure inverse data or idt on smoothed image?
    #use "faint" and "bright" distances together?
    
    inverse_data=np.where(distdata < 1.0e-3, np.ones_like(distdata), np.zeros_like(distdata))
    idt=scipy.ndimage.distance_transform_cdt(inverse_data)
    idt=scipy.ndimage.gaussian_filter(np.float32(idt),7)


    
    #print('200')
    #sm=np.where(idt > 200, scipy.ndimage.gaussian_filter(data,200), data)
    #print('100')
    #sm=np.where(idt > 100, scipy.ndimage.gaussian_filter(sm,100),data)
    #print('50')
    #sm=np.where(idt > 50.0, scipy.ndimage.gaussian_filter(data,100),data)

    print('25')
    sm=np.where(idt > 25.0, scipy.ndimage.gaussian_filter(data,50),data)

    
    print('10')
    sm=np.where(idt > 10.0, scipy.ndimage.gaussian_filter(sm,20),data)

    print('1')
    sm=np.where(idt > 1.0, scipy.ndimage.gaussian_filter(sm,8),data)
    sm=np.where(idt > 0.3, scipy.ndimage.gaussian_filter(sm,3),data)

    print('0.1')
    sm=np.where(idt > 0.1, scipy.ndimage.gaussian_filter(sm,2.5),data)
    sm=np.where(idt > 0.03, scipy.ndimage.gaussian_filter(sm,2),data)

    print('0')
    sm=np.where(idt > 0.0, scipy.ndimage.gaussian_filter(sm,1.5),data)

    sm=np.where(distdata < 1.0e-3, scipy.ndimage.gaussian_filter(data,2.5), sm)
    
    #do a final one for smoothness
    print('all')
    sm=scipy.ndimage.gaussian_filter(sm,1.0)
    
    outhdu=fits.PrimaryHDU(sm)
    outhdu.header=header

    disthdu=fits.ImageHDU(idt)
    
    outlist=fits.HDUList([outhdu,disthdu])
    outlist.writeto(new_fn,overwrite=True)
    
    hdulist.close()
    distlist.close()
    
    return


def adaptive_fast(sub_image,k_image,size_factor):

    inverse_data=np.where(k_image < (1.0e-3), np.ones_like(k_image), np.zeros_like(k_image))
    idt=scipy.ndimage.distance_transform_cdt(inverse_data)

    if sub_image.shape[0] > 200:
        idt=scipy.ndimage.gaussian_filter(np.float32(idt),20.0)
        
        sm=np.where(idt > 10.0, scipy.ndimage.gaussian_filter(sub_image,32.0),sub_image)                    
        sm=np.where(idt > 1.0, scipy.ndimage.gaussian_filter(sm,16.0),sub_image)
        sm=np.where(idt > 0.1, scipy.ndimage.gaussian_filter(sm,8.0),sub_image)
        sm=np.where(idt > 0.01, scipy.ndimage.gaussian_filter(sm,4.0),sub_image)
        sm=np.where(idt > 0.0, scipy.ndimage.gaussian_filter(sm,2.0),sub_image)
        sm=scipy.ndimage.gaussian_filter(sm,1.0)
                    
    elif sub_image.shape[0] > 100:
        idt=scipy.ndimage.gaussian_filter(np.float32(idt),10.0)
        
        sm=np.where(idt > 10.0, scipy.ndimage.gaussian_filter(sub_image,16.0),sub_image)        
        sm=np.where(idt > 1.0, scipy.ndimage.gaussian_filter(sm,8.0),sub_image)        
        sm=np.where(idt > 0.1, scipy.ndimage.gaussian_filter(sm,4.0),sub_image)        
        sm=np.where(idt > 0.01, scipy.ndimage.gaussian_filter(sm,2.0),sub_image)
        sm=scipy.ndimage.gaussian_filter(sm,1.0)


    elif sub_image.shape[0] > 50:
        idt=scipy.ndimage.gaussian_filter(np.float32(idt),5.0)
        
        sm=np.where(idt > 10.0, scipy.ndimage.gaussian_filter(sub_image,8.0),sub_image)        
        sm=np.where(idt > 1.0, scipy.ndimage.gaussian_filter(sm,4.0),sub_image)        
        sm=np.where(idt > 0.1, scipy.ndimage.gaussian_filter(sm,2.0),sub_image)        
        sm=scipy.ndimage.gaussian_filter(sm,1.0)

    elif sub_image.shape[0] > 25:
        idt=scipy.ndimage.gaussian_filter(np.float32(idt),3.0)

        sm=np.where(idt > 1.0, scipy.ndimage.gaussian_filter(sub_image,2.0),sub_image)        
        sm=np.where(idt > 0.1, scipy.ndimage.gaussian_filter(sm,1.0),sub_image)        
        sm=scipy.ndimage.gaussian_filter(sm,0.5)        

    elif sub_image.shape[0] > 10:
        idt=scipy.ndimage.gaussian_filter(np.float32(idt),2.0)

        sm=np.where(idt > 1.0, scipy.ndimage.gaussian_filter(sub_image,2.0),sub_image)        
        sm=np.where(idt > 0.1, scipy.ndimage.gaussian_filter(sm,1.0),sub_image)        
        sm=scipy.ndimage.gaussian_filter(sm,0.2)        

    else:
        sm=scipy.ndimage.gaussian_filter(sub_image,1.0)
        


    
     
    
    return sm

def adaptive_filter(values):

    l=values.shape[0]/2
    lrt=np.int32(l**0.5)
    
    reshape_im=np.reshape(values,(lrt,lrt,2))

    
    #im_bg=np.where(reshape_im[:,:,0]>1.0e-7,reshape_im[:,:,0],np.zeros_like(reshape_im[:,:,0]))

    #im_dt=scipy.ndimage.distance_transform_cdt(im_bg)

    im_dt=reshape_im[:,:,1]
    im_bg=reshape_im[:,:,0]
    
    half=np.int32(lrt/2)

    #distance=im_dt[half-1,half-1]

    #if distance > 3:
    #    returnval=im_bg[half-1,half-1]
    #else:
        #sigma_int=np.max( [ 5*(10-distance), 1] )
        #blur=scipy.ndimage.gaussian_filter(im_bg,sigma_int)
    #    returnval=np.mean(im_bg)

    #returnarr=np.where(im_dt < 10, np.mean(im_bg)*np.ones_like(im_bg), im_bg)

    #dt=float(im_dt[half-1,half-1])
    #weight=20.0/(1.0+dt)

    #val=im_bg[half-1,half-1]
    #expon=1.5
    #returnval= (np.mean(im_bg)*(weight)**expon + val)/(1.0+weight**expon)

    #too slow?
    #dval=np.int( 10.0*np.float( np.median(im_dt)+1 )**0.5 )
    #blur=scipy.ndimage.gaussian_filter(im_bg,dval)
    #returnval=blur[half-1,half-1]

    dval=np.int32(2.0*np.ceil( im_dt[half-1,half-1]) )  #already transformed

    i1=half-dval
    i2=half+dval
    i1=np.max([0,i1])
    i2=np.min([lrt-1,i2])

    if i1==i2:
        returnval=im_bg[i1,i1]
    else:
        returnval=np.mean(im_bg[i1:i2,i1:i2])  #way to do this with circular distance?
    #print(dval,half,lrt,i1,i2,returnval)
    
    return returnval

def process_single_filter_subimage(image_parameters,galaxy_data,lcdata,filname,lambda_eff_microns,selected_catalog=None):

    print('**** Processing subimage:  ',image_parameters)

    single_filter_subimage=np.ndarray((image_parameters['Npix'],image_parameters['Npix']))
    single_filter_subimage_smooth=np.ndarray((image_parameters['Npix'],image_parameters['Npix']))

    print(single_filter_subimage.shape)

    #find sources in this sub-image

    #use buffer to clear edge effects from lower left
    buf_deg=10.0/3600.0 #10 arcsec buffer?
    sub_indices=(lcdata['ra_deg']>=image_parameters['x1_deg']-buf_deg)*(lcdata['ra_deg']<image_parameters['x2_deg'])*(lcdata['dec_deg']>=image_parameters['y1_deg']-buf_deg)*(lcdata['dec_deg']<image_parameters['y2_deg'])

    sub_data=lcdata[sub_indices]

    success=0

    #image_catalog={'filter':filname}

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


    number=len(sub_data)

    if selected_catalog is not None:
        assert(len(selected_catalog['galaxy_indices'])==number)
        build_catalog=False
    else:
        selected_catalog={'filter':filname}
        selected_catalog['galaxy_indices']=np.ndarray((number),dtype=object)
        selected_catalog['found_galaxy']=np.ndarray((number),dtype=object)
        selected_catalog['this_camstr']=np.ndarray((number),dtype=object)
        selected_catalog['simlabel']=np.ndarray((number),dtype=object)
        selected_catalog['numcams']=np.ndarray((number),dtype=object)
        selected_catalog['image_dir']=np.ndarray((number),dtype=object)
        selected_catalog['scalefactor']=np.ndarray((number),dtype=object)
        selected_catalog['Mvir']=np.ndarray((number),dtype=object)
        selected_catalog['mstar']=np.ndarray((number),dtype=object)
        selected_catalog['rhalf']=np.ndarray((number),dtype=object)
        selected_catalog['lc_entry']=np.ndarray((number),dtype=object)
        selected_catalog['orig_pix_arcsec']=np.ndarray((number),dtype=object)

        selected_catalog['mstar_factor']=np.ndarray((number),dtype=object)
        selected_catalog['size_factor']=np.ndarray((number),dtype=object)
        
        selected_catalog['icen']=np.ndarray((number),dtype=object)
        selected_catalog['jcen']=np.ndarray((number),dtype=object)
        selected_catalog['flux_njy']=np.ndarray((number),dtype=object)
        selected_catalog['ABmag']=np.ndarray((number),dtype=object)
        selected_catalog['in_image']=np.ndarray((number),dtype=object)

        selected_catalog['final_file']=np.ndarray((number),dtype=object)
        selected_catalog['Kband_file']=np.ndarray((number),dtype=object)
        

        
        build_catalog=True

    
    data_found=0

    for i,entry in enumerate(sub_data):

        if i % 100 == 0:
            print('   processed ', i, ' out of ', number)
            
        #need pos_i, pos_j
        pos_i=np.int64( (entry['ra_deg']-image_parameters['x1_deg'])*np.float64(image_parameters['Npix'])/(image_parameters['x2_deg']-image_parameters['x1_deg'])   )
        pos_j=np.int64( (entry['dec_deg']-image_parameters['y1_deg'])*np.float64(image_parameters['Npix'])/(image_parameters['y2_deg']-image_parameters['y1_deg'])   )

        #select image file to insert
        mass_value=entry['subhalo_mass_msun']
        scale_value=1.0/(1.0+entry['true_z'])
        mstar_value=entry['mstar_msun_rad']
        

        if build_catalog is True:
            galaxy_scale_indices=(galaxy_data['scale']>=scale_value/scale_window)
            galaxy_scale_indices*=(galaxy_data['scale']<=scale_value*scale_window)
            
            galaxy_mass_indices=(galaxy_data['Mvir']>=mass_value/mass_window)
            galaxy_mass_indices*=(galaxy_data['Mvir']<=mass_value*mass_window)
            
            galaxy_search_indices=galaxy_scale_indices*galaxy_mass_indices
            
            found_galaxy=False

            selected_catalog['lc_entry'][i]=entry

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

            selected_catalog['found_galaxy'][i]=found_galaxy
            selected_catalog['galaxy_indices'][i]=galaxy_index
        else:
            galaxy_index=selected_catalog['galaxy_indices'][i]
            found_galaxy=selected_catalog['found_galaxy'][i]


        
        found_data=False
        #now we have galaxy_index
        if galaxy_index is not None:
            mstar_factor=mstar_value/(galaxy_data['Mstar'][galaxy_index])
            size_factor=(mstar_factor)**0.5

            folder=galaxy_data['image_dir'][galaxy_index]
            label=galaxy_data['simlabel'][galaxy_index]
            if build_catalog is True:
                selected_catalog['simlabel'][i]=label
                selected_catalog['image_dir'][i]=folder

                possible_files=np.sort(np.asarray(glob.glob(folder+'/hires_images_cam??/'+label+'cam??_'+filname+'*.fits')))
                selected_catalog['numcams'][i]=possible_files.shape[0]
                
                if possible_files.shape[0]>0:
                    #pick a random camera
                    file_index=random.randint(possible_files.shape[0])
                    #assert all filters exist
                    this_file=possible_files[file_index]
                    this_folder=os.path.dirname(this_file)
                    this_camstr=this_folder[-5:]
                    filter_files=np.sort(np.asarray(glob.glob(this_folder+'/'+label+'*.fits')))
                    if filter_files.shape[0]==8:
                        kband_files=np.asarray(glob.glob(folder+'/hires_images_'+this_camstr+'/'+label+this_camstr+'_'+'nircam_f200w*.fits'))
                        if kband_files.shape[0]==0:
                            assert(false)
                        else:
                            selected_catalog['Kband_file'][i]=kband_files[0]
                    else:
                        this_camstr=None
                        
                    selected_catalog['this_camstr'][i]=this_camstr
                else:
                    this_file=None
            else:
                this_camstr=selected_catalog['this_camstr'][i]
                if this_camstr is not None:
                    this_files=np.asarray(glob.glob(folder+'/hires_images_'+this_camstr+'/'+label+this_camstr+'_'+filname+'*.fits'))
                    if this_files.shape[0]==0:
                        assert(False)
                    this_file=this_files[0]

                else:
                    this_file=None


            selected_catalog['final_file'][i]=this_file
            
            if this_file is not None:
                found_data=True
                
                this_hdu=fits.open(this_file)[0]
                this_image=this_hdu.data
                
                pixstr=this_file.split('_')[-2][3:]
                
                pixsize_arcsec=np.float64(pixstr)

                #adaptive-smmoother here???  HOW?
                kband_hdu=fits.open(selected_catalog['Kband_file'][i])[0]
                kdata=kband_hdu.data

                #measure inverse? distance transform in K filter
                #not sure this is the best flux limit.. a higher one will give more agressive blurring in outer gal
                
                #inverse_im=np.where(kdata < 1.0e-4, np.ones_like(kdata), np.zeros_like(kdata))
                #idt=scipy.ndimage.distance_transform_cdt(inverse_im)

                #run generic_filter with adaptive_filter function
                #inarr=np.ndarray((kdata.shape[0],kdata.shape[0],2),dtype=np.float64)
                #outarr=np.zeros_like(inarr)
                #inarr[:,:,0]=kdata
                #inarr[:,:,1]=idt  #transform this

                #actually.. apply to re-sized images for efficiency???  need ultra-fast image transforms..
                
            else:
                found_data=False


            if found_data==True:
                #these are in nJy-- preserve integral!
                data_found+=1
                
                original_flux=np.sum(this_image)
                
                total_flux=original_flux*mstar_factor  #shrink luminosity by mstar factor
                
                this_npix=this_image.shape[0]
                
                #resize--preserve proper units
                
                desired_npix=np.int32( this_npix*(pixsize_arcsec/image_parameters['pix_size_arcsec'])*size_factor  )
                
                resized_image=congrid.congrid(this_image,(desired_npix,desired_npix))
                resized_flux=np.sum(resized_image)
                
                resized_image=resized_image*(total_flux/resized_flux)
                
                resized_k=congrid.congrid(kdata,(desired_npix,desired_npix))

                #inverse_k=np.where(resized_k < 1.0e-3,np.ones_like(resized_k),np.zeros_like(resized_k))
                #idt=scipy.ndimage.distance_transform_cdt(inverse_k)
                #inarr=np.ndarray((idt.shape[0],idt.shape[0],2),dtype=np.float64)
                #outarr=np.zeros_like(inarr)

                #inarr[:,:,0]=resized_image
                #inarr[:,:,1]=scipy.ndimage.gaussian_filter(2.5*idt**0.5,5)

                #scipy.ndimage.generic_filter(inarr,adaptive_filter,size=(10,10,2),output=outarr,origin=0.5)

                #tform_image=outarr[:,:,1]
                
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


                #I think this is wrong?
                #sub_image_to_add=resized_image[im0-i1:npsub-(i2-im1),im0-j1:npsub-(j2-im1)]
                #k_subimage=resized_k[im0-i1:npsub-(i2-im1),im0-j1:npsub-(j2-im1)]


                if i1 < im0:
                    is1=np.abs(im0-i1)
                else:
                    is1=0

                if i2 >= im1:
                    is2=npsub-(i2-im1)-1
                else:
                    is2=npsub-1

                if j1 < im0:
                    js1=np.abs(im0-j1)
                else:
                    js1=0

                if j2 >= im1:
                    js2=npsub-(j2-im1)-1
                else:
                    js2=npsub-1


                sub_image_to_add=resized_image[is1:is2,js1:js2]
                k_subimage=resized_k[is1:is2,js1:js2]
                    
                new_image_to_add=adaptive_fast(sub_image_to_add,k_subimage,size_factor)
                k_smoothed=adaptive_fast(k_subimage,k_subimage,size_factor)
                orig_image_to_add=sub_image_to_add
                
                #detect edges?

                if k_subimage.shape[0] > 2 and k_subimage.shape[1] > 2:
                    edge1=np.mean(k_smoothed[0:2,:])
                    edge2=np.mean(k_smoothed[:,-2:])
                    edge3=np.mean(k_smoothed[:,0:2])
                    edge4=np.mean(k_smoothed[-2:,:])
                    maxedge=np.max(np.asarray([edge1,edge2,edge3,edge4]))
                else:
                    maxedge=0.0


                #print('edge ratios: ', edge1/total_flux, edge2/total_flux, edge3/total_flux, edge4/total_flux, ' file ', this_file, size_factor)

                in_image=False
                
                if maxedge > 0.03 and resized_image.shape[0] > 20 and new_image_to_add.shape[0]==new_image_to_add.shape[1]:
                    print('omitting edge effect, max: ', maxedge, os.path.basename(this_file), size_factor, new_image_to_add.shape)
                    new_image_to_add *= 0.0
                    orig_image_to_add *= 0.0
                    in_image=False
                elif maxedge > 0.001 and resized_image.shape[0] > 200 and new_image_to_add.shape[0]==new_image_to_add.shape[1]:
                    print('omitting edge effect, max: ', maxedge, os.path.basename(this_file), size_factor, new_image_to_add.shape)
                    new_image_to_add *= 0.0
                    orig_image_to_add *= 0.0
                    in_image=False                    
                else:
                    in_image=True
                    
                if icen >= 0.0 and jcen >= 0.0:
                    #assemble catalog entries for this object
                    #final source info
                    xcen_list.append(icen)
                    ycen_list.append(jcen)
                    
                    final_flux_njy_list.append(np.sum(new_image_to_add))
                    if np.sum(new_image_to_add)==0.0:
                        ab_appmag_list.append(-1)
                    else:    
                        ab_appmag_list.append(-2.5*np.log10((1.0e9)*(np.sum(new_image_to_add))/3632.0))
                else:
                    in_image=False


                    
                if sub_image_to_add.shape[0] > 0 and sub_image_to_add.shape[1] > 0:
                    iplace=np.max([i1,0])
                    jplace=np.max([j1,0])
                    new_npi=sub_image_to_add.shape[0]
                    new_npj=sub_image_to_add.shape[1]
                    single_filter_subimage[iplace:iplace+new_npi,jplace:jplace+new_npj] += orig_image_to_add
                    single_filter_subimage_smooth[iplace:iplace+new_npi,jplace:jplace+new_npj] += new_image_to_add
                    


                        
                selected_catalog['scalefactor'][i]=galaxy_data['scale'][galaxy_index]
                selected_catalog['Mvir'][i]=galaxy_data['Mvir'][galaxy_index]
                selected_catalog['mstar'][i]=galaxy_data['Mstar'][galaxy_index]
                selected_catalog['rhalf'][i]=galaxy_data['Rhalf_stars'][galaxy_index]
                selected_catalog['orig_pix_arcsec'][i]=pixsize_arcsec

                selected_catalog['mstar_factor'][i]=mstar_factor
                selected_catalog['size_factor'][i]=size_factor

                
                selected_catalog['icen'][i]=icen
                selected_catalog['jcen'][i]=jcen
                selected_catalog['flux_njy'][i]=total_flux
                selected_catalog['ABmag'][i]=-2.5*np.log10((1.0e9)*(total_flux)/3632.0)
                selected_catalog['in_image'][i]=in_image


            
            
        
    print('**** Subimage data found: ', str(data_found), '  out of  ', i+1)

    
    return single_filter_subimage,single_filter_subimage_smooth,selected_catalog


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
    
    filters = np.flipud(['U','B','V','Z','nircam_f115w','nircam_f150w','nircam_f200w'])
    output_filters=np.flipud(['f336w','f435w','f606w','f850lp','f115w','f150w','f200w'])
    lambda_eff_microns = np.flipud([0.35,0.45,0.55,0.85,1.15,1.50,2.0])


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

    selected_catalog=None

    #pick random seed
    baseseed=10

    random.seed(baseseed+image_parameters['xsub']+int(total_images**0.5)*image_parameters['ysub'])

    for i,filname in enumerate(filters):
        print('processing.. ', filname)


        
        single_filter_subimage,single_filter_subimage_smooth,selected_catalog=process_single_filter_subimage(image_parameters,
                                                                                             galaxy_data,
                                                                                             lcdata,
                                                                                             filname,
                                                                                             lambda_eff_microns[i],
                                                                                             selected_catalog=selected_catalog)

        output_filename=os.path.join(output_dir,image_filelabel+'_imager'+'_fielda-subimage'+str(do_sub_image_x)+str(do_sub_image_y)+'-9-8_'+output_filters[i]+'_v1_lightcone.fits' )
        output_filename_sm=os.path.join(output_dir,image_filelabel+'_imager'+'_fielda-subimage'+str(do_sub_image_x)+str(do_sub_image_y)+'-9-8-smooth_'+output_filters[i]+'_v1_lightcone.fits' )

        print('**** Saving Subimage: ', output_filename )

        #print(selected_catalog)

        #save this catalog somehow???
        
        outhdu=fits.PrimaryHDU(np.float32(single_filter_subimage))
        outhdu.header['PIXSIZE']=(image_parameters['pix_size_arcsec'],'arcsec')
        
        outlist=fits.HDUList([outhdu])
        outlist.writeto(output_filename,overwrite=True)


        outhdu2=fits.PrimaryHDU(np.float32(single_filter_subimage_smooth))
        outhdu2.header['PIXSIZE']=(image_parameters['pix_size_arcsec'],'arcsec')
        
        outlist2=fits.HDUList([outhdu2])
        outlist2.writeto(output_filename_sm,overwrite=True)
        
        
        
    return
