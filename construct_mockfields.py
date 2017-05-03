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

sq_arcsec_per_sr = 42545170296.0
c = 3.0e8

filters_to_analyze = np.asarray(['hst/acs_f435w','hst/acs_f606w','hst/acs_f775w','hst/acs_f850lp',
                      'hst/wfc3_f105w','hst/wfc3_f125w','hst/wfc3_f160w',
                      'jwst/nircam_f070w', 'jwst/nircam_f090w','jwst/nircam_f115w', 'jwst/nircam_f150w', 
                      'jwst/nircam_f200w', 'jwst/nircam_f277w', 'jwst/nircam_f356w', 'jwst/nircam_f444w', 
                      'hst/wfc3_f140w',
                      'hst/wfc3_f275w', 'hst/wfc3_f336w',
                      'hst/acs_f814w',
                      'jwst/miri_F560W','jwst/miri_F770W','jwst/miri_F1000W','jwst/miri_F1130W',
                      'jwst/miri_F1280W','jwst/miri_F1500W','jwst/miri_F1800W','jwst/miri_F2100W','jwst/miri_F2550W'])

psf_dir = os.path.expandvars('$GFS_PYTHON_CODE/vela-yt-sunrise/kernels')

psf_names = np.asarray(['TinyTim_IllustrisPSFs/F435W_rebin.fits','TinyTim_IllustrisPSFs/F606W_rebin.fits','TinyTim_IllustrisPSFs/F775W_rebin.fits','TinyTim_IllustrisPSFs/F850LP_rebin.fits',
             'TinyTim_IllustrisPSFs/F105W_rebin.fits','TinyTim_IllustrisPSFs/F125W_rebin.fits','TinyTim_IllustrisPSFs/F160W_rebin.fits',
             'WebbPSF_F070W_trunc.fits','WebbPSF_F090W_trunc.fits','WebbPSF_F115W_trunc.fits','WebbPSF_F150W_trunc.fits',
             'WebbPSF_F200W_trunc.fits','WebbPSF_F277W_trunc.fits','WebbPSF_F356W_trunc.fits','WebbPSF_F444W_trunc.fits',
             'TinyTim_IllustrisPSFs/F140W_rebin.fits','TinyTim_IllustrisPSFs/F275W_rebin.fits','TinyTim_IllustrisPSFs/F336W_rebin.fits','TinyTim_IllustrisPSFs/F814W_rebin.fits',
             'WebbPSF_F560W_trunc.fits','WebbPSF_F770W_trunc.fits','WebbPSF_F1000W_trunc.fits','WebbPSF_F1130W_trunc.fits',
             'WebbPSF_F1280W_trunc.fits','WebbPSF_F1500W_trunc.fits','WebbPSF_F1800W_trunc.fits','WebbPSF_F2100W_trunc.fits','WebbPSF_F2550W_trunc.fits'])


psf_pix_arcsec = np.asarray([0.03,0.03,0.03,0.03,0.06,0.06,0.06,0.0317,0.0317,0.0317,0.0317,0.0317,0.0648,0.0648,0.0648,0.06,0.03,0.03,0.03,0.11,0.11,0.11,0.11,0.11,0.11,0.11,0.11,0.11])
psf_fwhm = np.asarray([0.10,0.11,0.12,0.13,0.14,0.17,0.20,0.11,0.11,0.11,0.11,0.12,0.15,0.18,0.25,0.18,0.07,0.08,0.13,
            0.035*5.61,0.035*7.57,0.035*9.90,0.035*11.30,0.035*12.75,0.035*14.96,0.035*17.90,0.035*20.65,0.035*25.11])

#photfnu units Jy; flux in 1 ct/s
photfnu_Jy = np.asarray([1.96e-7,9.17e-8,1.97e-7,4.14e-7,
                         1.13e-7,1.17e-7,1.52e-7,
                         5.09e-8,3.72e-8,3.17e-8,2.68e-8,2.64e-8,2.25e-8,2.57e-8,2.55e-8,
                         9.52e-8,8.08e-7,4.93e-7,1.52e-7,
                         5.75e-8,3.10e-8,4.21e-8,1.39e-7,
                         4.65e-8,4.48e-8,5.88e-8,4.98e-8,1.15e-7])

#construct real illustris lightcones from individual images

#in parallel, produce estimated Hydro-ART surveys based on matching algorithms -- high-res?


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


def process_single_filter(data,lcdata,filname,fil_index,output_dir,image_filelabel,image_suffix,eff_lambda_microns,lim=None,minz=None):

    data=copy.copy(data)
    
    print('Processing:  ', filname)
    full_npix=data['full_npix'][0]
    pixsize_arcsec=data['pixsize_arcsec'][0]
    n_galaxies=data['full_npix'].shape[0]


    if filname.find('WFI')==0:
        fn=filname[-4:]
        filname='wfirst/wfidrm15_'+fn
        
    try:
        pbi= filters_to_analyze==filname
        this_psf_file=os.path.join(psf_dir,psf_names[pbi][0])
        this_psf_pixsize_arcsec=psf_pix_arcsec[pbi][0]
        this_psf_fwhm=psf_fwhm[pbi][0]
        this_photfnu_Jy=photfnu_Jy[pbi][0]
        print('PSF info: ', this_psf_file, this_psf_pixsize_arcsec, this_psf_fwhm, this_photfnu_Jy)
        do_psf=True
    except:
        print('Missing filter info, skipping PSF: ', filname)
        do_psf=False
        this_psf_pixsize_arcsec=pixsize_arcsec
        #return None
    


    desired_pixsize_arcsec=this_psf_pixsize_arcsec

    full_fov=full_npix*pixsize_arcsec

    desired_npix=full_fov/desired_pixsize_arcsec

    print('Orig pix: ', full_npix, ' Desired pix: ', desired_npix)

    if do_psf is True:
        orig_psf_kernel = pyfits.open(this_psf_file)[0].data 

        #psf kernel shape must be odd for astropy.convolve??
        if orig_psf_kernel.shape[0] % 2 == 0:
            new_psf_shape = orig_psf_kernel.shape[0]-1
            psf_kernel = congrid.congrid(orig_psf_kernel,(new_psf_shape,new_psf_shape))
        else:
            psf_kernel = orig_psf_kernel
        
        assert( psf_kernel.shape[0] % 2 != 0)


    image_cube = np.zeros((full_npix,full_npix),dtype=np.float64)

    success=[]
    mag=[]

    #for bigger files, may need to split by filter first
    index=np.arange(n_galaxies)

    for pos_i,pos_j,origin_i,origin_j,run_dir,this_npix,this_z,num in zip(data['pos_i'],data['pos_j'],data['origin_i'],data['origin_j'],data['run_dir'],data['this_npix'],data['z'],index):
        if lim is not None:
            if num > lim:
                success.append(False)
                continue
        if minz is not None:
            if this_z < minz:
                success.append(False)
                continue

        try:
            bblist=pyfits.open(os.path.join(run_dir,'broadbandz.fits'))
            this_cube = bblist['CAMERA0-BROADBAND-NONSCATTER'].data
            this_mag = (bblist['FILTERS'].data['AB_mag_nonscatter0'])[fil_index]
            bblist.close()

            #if catalog and image have different npix, this is a failure somewhere
            cube_npix=this_cube.shape[-1]
            assert(cube_npix==this_npix)
            
            success.append(True)
            mag.append(this_mag)
        except:
            print('Missing file or mismatched shape, ', run_dir, cube_npix, this_npix)
            success.append(False)
            continue



            
        i_tc=0
        j_tc=0
        i_tc1=this_npix
        j_tc1=this_npix

        if origin_i < 0:
            i0=0
            i_tc=-1*origin_i
        else:
            i0=origin_i

        if origin_j < 0:
            j0=0
            j_tc=-1*origin_j
        else:
            j0=origin_j

        if i0+this_npix > full_npix:
            i1=full_npix
            i_tc1= full_npix-i0   #this_npix - (i0+this_npix-full_npix)
        else:
            i1=i0+this_npix-i_tc

        if j0+this_npix > full_npix:
            j1=full_npix
            j_tc1= full_npix-j0
        else:
            j1=j0+this_npix-j_tc


        sub_cube1=image_cube[i0:i1,j0:j1]
        this_subcube=this_cube[fil_index,i_tc:i_tc1,j_tc:j_tc1]
        print(run_dir, this_subcube.shape)

        image_cube[i0:i1,j0:j1] = sub_cube1 + this_subcube



    #convolve here
    
    #first, re-grid to desired scale

    new_image=congrid.congrid(image_cube,(desired_npix,desired_npix))
    new_i=data['pos_i']*desired_npix/full_npix
    new_j=data['pos_j']*desired_npix/full_npix

    
    pixel_Sr = (desired_pixsize_arcsec**2)/sq_arcsec_per_sr  #pixel area in steradians:  Sr/pixel
    to_nJy_per_Sr = (1.0e9)*(1.0e14)*(eff_lambda_microns**2)/c   #((pixscale/206265.0)^2)*
    #sigma_nJy = 0.3*(2.0**(-0.5))*((1.0e9)*(3631.0/5.0)*10.0**(-0.4*self.maglim))*self.Pix_arcsec*(3.0*self.FWHM_arcsec)
    to_nJy_per_pix = to_nJy_per_Sr*pixel_Sr
        
    nopsf_im=new_image*to_nJy_per_pix

    if do_psf is True:
        conv_im = convolve_fft(new_image,psf_kernel,boundary='fill',fill_value=0.0,normalize_kernel=True,allow_huge=True)
        final_im=conv_im*to_nJy_per_pix
    
    outname=os.path.join(output_dir,image_filelabel+'_'+filname.replace('/','-')+'_'+image_suffix+'_v1_lightcone.fits')
    print('saving:', outname)

    primary_hdu=pyfits.PrimaryHDU(nopsf_im)
    primary_hdu.header['FILTER']=filname.replace('/','-')
    primary_hdu.header['PIXSIZE']=(desired_pixsize_arcsec,'arcsec')
    primary_hdu.header['UNIT']=('nanoJanskies','per pixel')
    abzp= - 2.5*(-9.0) + 2.5*np.log10(3631.0)  #images in nanoJanskies
    primary_hdu.header['ABZP']=(abzp, 'AB mag zeropoint')
    if do_psf is True:
        primary_hdu.header['PHOTFNU']=(this_photfnu_Jy,'Jy; approx flux[Jy] at 1 count/sec')
        
    primary_hdu.header['EXTNAME']='IMAGE_NOPSF'

    if do_psf is True:
        psfim_hdu=pyfits.ImageHDU(final_im)
        psfim_hdu.header['EXTNAME']='IMAGE_PSF'
        
        psf_hdu = pyfits.ImageHDU(psf_kernel)
        psf_hdu.header['EXTNAME']='MODELPSF'
        psf_hdu.header['PIXSIZE']=(desired_pixsize_arcsec,'arcsec')


    if np.sum(np.asarray(data.colnames)=='success')==0:
        newcol=astropy.table.column.Column(data=success,name='success')
        data.add_column(newcol)

    if np.sum(np.asarray(data.colnames)=='new_i')==0:
        newicol=astropy.table.column.Column(data=new_i,name='new_i')
        newjcol=astropy.table.column.Column(data=new_j,name='new_j')

        data.add_column(newicol)
        data.add_column(newjcol)

    magcol=astropy.table.column.Column(data=mag,name='AB_absmag_'+filname.replace('/','-'))
    data.add_column(magcol)
    
    data_df=data.to_pandas()
    lc_df = lcdata.to_pandas()
    lc_df.rename(columns=lcfile_cols,inplace=True)
    
    assert(lc_df.shape[0]==data_df.shape[0])
    
    new_df=lc_df.join(data_df)

    failures = new_df.where(new_df['success']==False).dropna()
    successes=new_df.drop(failures.index)
    print('N successes: ', successes.shape[0])
    
    new_data=astropy.table.Table.from_pandas(successes)



    table_hdu = pyfits.table_to_hdu(new_data)
    table_hdu.header['EXTNAME']='Catalog'

    if do_psf is True:
        output_list=pyfits.HDUList([primary_hdu,psfim_hdu,psf_hdu,table_hdu])
    else:
        output_list=pyfits.HDUList([primary_hdu,table_hdu])



    tempfile=os.path.join(os.path.expandvars('/scratch/$USER/$SLURM_JOBID'),os.path.basename(outname))
    print('saving to scratch first.. , ', tempfile)
    output_list.writeto(tempfile,overwrite=True)
    output_list.close()

    shutil.copy(tempfile,output_dir)
    
    return success


def build_lightcone_images(image_info_file,lightcone_file,run_type='images',lim=None,minz=None,image_filelabel='hlsp_misty_illustris',jwst_only=False):

    data=ascii.read(image_info_file)
    print(data)


    #get expected shape
    test_file=os.path.join(data['run_dir'][0],'broadbandz.fits')
    tfo =pyfits.open(test_file)
    print(tfo.info())

    try:
        cube=tfo['CAMERA0-BROADBAND-NONSCATTER'].data
        cubeshape=cube.shape
    except:
        square=tfo['CAMERA0-PARAMETERS'].data
        nf=tfo['FILTERS'].data.shape[0]
        cubeshape=(nf,square.shape[0],square.shape[0])
        
    print(cubeshape)
    auxcube=tfo['CAMERA0-AUX'].data

    filters_hdu = tfo['FILTERS']

    lightcone_dir=os.path.abspath(os.path.dirname(image_info_file))
    print('Constructing lightcone data from: ', lightcone_dir)

    output_dir = os.path.join(lightcone_dir,os.path.basename(image_info_file).rstrip('.txt'))
    print('Saving lightcone outputs in: ', output_dir)
    if not os.path.lexists(output_dir):
        os.mkdir(output_dir)


    image_suffix=os.path.basename(image_info_file).rstrip('_images.txt')
        
    success_catalog=os.path.join(output_dir,os.path.basename(image_info_file).rstrip('.txt')+'_success.txt')

    N_filters = cubeshape[0]

    #N_aux=auxcube.shape[0]
    #aux_cube = np.zeros((N_aux,full_npix,full_npix),dtype=np.float64)

    lcdata=ascii.read(lightcone_file)

    filters_data=filters_hdu.data
    print(filters_data.columns)
    lambda_eff_microns = filters_data['lambda_eff']*1.0e6

    for i,filname in enumerate(filters_data['filter']):
        if jwst_only is True:
            if filname.find('jwst')==-1:
                continue

        success=process_single_filter(data,lcdata,filname,i,output_dir,image_filelabel,image_suffix,lambda_eff_microns[i],lim=lim,minz=minz)
        if i==0:
            success=np.asarray(success)
            #newcol=astropy.table.column.Column(data=success,name='success')
            #data.add_column(newcol)
            ascii.write(data,output=success_catalog,overwrite=True)


    #convert units before saving.. or save both?

    return
