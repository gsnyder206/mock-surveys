import astropy
import astropy.io.fits as fits
import numpy as np
import gfs_sublink_utils as gsu
import make_color_image
import matplotlib
import matplotlib.pyplot as pyplot
import glob
import os
import time
import sys

sq_arcsec_per_sr = 42545170296.0
c = 3.0e8


def export_image(hdulist,camnum,filtername,label='',outdir='.',nonscatter=False,filterlabel=None):

    #get image in W/m/m^2/Sr

    #hdulist=fits.open(bbfile)
    fils=hdulist['FILTERS'].data['filter']
    fi=np.where(fils==filtername)[0][0]
    
    efl=hdulist['FILTERS'].data['lambda_eff']
    
    efl_microns=1.0e6 * efl[fi]
    
    if nonscatter is False:
        key='CAMERA'+'{}'.format(camnum)+'-BROADBAND'
    else:
        key='CAMERA'+'{}'.format(camnum)+'-BROADBAND-NONSCATTER'

        

    print('loading data into variable.. ')
    s=time.time()
    camdata=hdulist[key].data[fi,:,:]
    f=time.time()
    print('loading data finished in (s): ',f-s )
    
    print(camdata.shape)
    #convert to nJy

    redshift=hdulist['BROADBAND'].header['redshift']
    pixsize_kpc=hdulist[key].header['CD1_1']
    pixsize_arcsec=pixsize_kpc/(gsu.illcos.kpc_proper_per_arcmin(redshift).value/60.0)
    
    print(pixsize_kpc,pixsize_arcsec)
    
    pixel_Sr = (pixsize_arcsec**2)/sq_arcsec_per_sr  #pixel area in steradians:  Sr/pixel
    to_nJy_per_Sr = (1.0e9)*(1.0e14)*(efl_microns**2)/c   #((pixscale/206265.0)^2)*
    to_nJy_per_pix = to_nJy_per_Sr*pixel_Sr

    newdata=camdata*to_nJy_per_pix

    
    #save new image
    if filterlabel is None:
        fname=filtername.split('/')[1]
    else:
        fname=filterlabel
        
    outputname=os.path.join(outdir, label+'cam'+'{:02d}'.format(camnum)+'_'+fname+'_pix'+'{:6.4f}'.format(pixsize_arcsec)+'_nJy.fits')
    print(outputname)

    outhdu=fits.PrimaryHDU(newdata)
    outhdu.header['redshift']=hdulist['BROADBAND'].header['redshift']
    outhdu.header['PIXSIZE']=(pixsize_arcsec,'arcsec')
    outhdu.header['PIXKPC']=(pixsize_kpc,'kpc')
    
    outlist=fits.HDUList([outhdu])
    outlist.writeto(outputname,overwrite=True)
    
    #fits.writeto(outputname,newdata,overwrite=True)
    
    return newdata,outputname


def export_hdst_filters(hdu_obj,camnum,vl,outdir='.',nonscatter=False):

    if hdu_obj.__class__ != fits.hdu.hdulist.HDUList:
        #pass filenames instead of hdu obj
        these_files=hdu_obj
        fn_u=these_files[np.where(np.core.defchararray.find(these_files, 'cam{:02d}'.format(camnum)+'_U_') != -1)[0]][0] ; u=(fits.open(fn_u)[0]).data
        fn_b=these_files[np.where(np.core.defchararray.find(these_files, 'cam{:02d}'.format(camnum)+'_B_') != -1)[0]][0] ; b=fits.open(fn_b)[0].data
        fn_v=these_files[np.where(np.core.defchararray.find(these_files, 'cam{:02d}'.format(camnum)+'_V_') != -1)[0]][0] ; v=fits.open(fn_v)[0].data
        fn_i=these_files[np.where(np.core.defchararray.find(these_files, 'cam{:02d}'.format(camnum)+'_I_') != -1)[0]][0] ; i=fits.open(fn_i)[0].data
        fn_z=these_files[np.where(np.core.defchararray.find(these_files, 'cam{:02d}'.format(camnum)+'_Z_') != -1)[0]][0] ; z=fits.open(fn_z)[0].data

        fn_j=these_files[np.where(np.core.defchararray.find(these_files, 'cam{:02d}'.format(camnum)+'_nircam_f115w_') != -1)[0]][0] ; j=fits.open(fn_j)[0].data
        fn_h=these_files[np.where(np.core.defchararray.find(these_files, 'cam{:02d}'.format(camnum)+'_nircam_f150w_') != -1)[0]][0] ; h=fits.open(fn_h)[0].data
        fn_k=these_files[np.where(np.core.defchararray.find(these_files, 'cam{:02d}'.format(camnum)+'_nircam_f200w_') != -1)[0]][0] ; k=fits.open(fn_k)[0].data

        redshift=fits.open(fn_v)[0].header['redshift']
        
    else:
        u,fn_u=export_image(hdu_obj,camnum,'hst/wfc3_f336w',label=vl,outdir=outdir,filterlabel='U')
        b,fn_b=export_image(hdu_obj,camnum,'hst/acs_f435w',label=vl,outdir=outdir,filterlabel='B')
        v,fn_v=export_image(hdu_obj,camnum,'hst/acs_f606w',label=vl,outdir=outdir,filterlabel='V')
        i,fn_i=export_image(hdu_obj,camnum,'hst/acs_f775w',label=vl,outdir=outdir,filterlabel='I')
        z,fn_z=export_image(hdu_obj,camnum,'hst/acs_f850lp',label=vl,outdir=outdir,filterlabel='Z')
        
        j,fn_j=export_image(hdu_obj,camnum,'jwst/nircam_f115w',label=vl,outdir=outdir,filterlabel='nircam_f115w')
        h,fn_h=export_image(hdu_obj,camnum,'jwst/nircam_f150w',label=vl,outdir=outdir,filterlabel='nircam_f150w')
        k,fn_k=export_image(hdu_obj,camnum,'jwst/nircam_f200w',label=vl,outdir=outdir,filterlabel='nircam_f200w')

        redshift=hdu_obj['BROADBAND'].header['redshift']


    scalefactor=(1.0/(1.0+redshift))
    
    rr=h+k
    gg=i+z+j
    bb=u+b+v

    fig = pyplot.figure(figsize=(6,6), dpi=600)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)
  
    axi = fig.add_subplot(1,1,1)
    axi.set_xticks([]) ; axi.set_yticks([])
    
    alph=10.0*(0.33/scalefactor)**5 ; Q=5.0
    
    rgbthing = make_color_image.make_interactive_nasa(bb,gg,rr,alph,Q)
    axi.imshow(rgbthing,interpolation='nearest',aspect='auto',origin='lower')

    outfile=os.path.join(outdir,'rgb_'+vl+'.pdf')
    
    fig.savefig(outfile,dpi=300,facecolor='Black')
    pyplot.close(fig)
    
    return


def export_run(sim_folder,camnums=np.arange(19),nonscatter=False,outdir='/astro/snyder_lab2/New_HydroART_images/VELA_v2/luvoir_mocks'):

    sunrise_dirs=np.sort(np.asarray(glob.glob(os.path.join(sim_folder,'*_sunrise'))))
    print(sunrise_dirs)

    sim_label=os.path.basename(sim_folder)
    print(sim_label)

    output_dir=os.path.join(outdir,sim_label)
    print(output_dir)

    if not os.path.lexists(output_dir):
        os.mkdir(output_dir)

    

    for sim_dir in sunrise_dirs:
        time_folder=os.path.basename(sim_dir)
        time_label=time_folder.split('_')[-2]
        print(time_label)

        bb_file=os.path.join(sim_dir,'images/broadbandz.fits.gz')
        if not os.path.lexists(bb_file):
            bb_file=os.path.join(sim_dir,'images/broadbandz.fits')
            if not os.path.lexists(bb_file):
                continue
            
        #check if exported fits files exist already.. in export_hdst_filters_function ?
        these_files=np.asarray(glob.glob(os.path.join(output_dir,'*/'+sim_label+'_'+time_label+'*.fits')),dtype=np.str_)

        print(these_files.shape, these_files, bb_file, sim_dir, output_dir, sim_label, time_label)

        if these_files.shape[0]==8*camnums.shape[0]:
            print('Files exist: ', these_files.shape[0], these_files[0])
            hdu_obj=these_files
        else:
            #print('loading hdulist... ')
            #s=time.time()

            hdu_obj=fits.open(bb_file,lazy_load_hdus=False)

            #load all camera data 

            
            #this doesn't work:  hdu_obj[slice:slice].readall()
            l=len(hdu_obj)
            
            #f=time.time()
            #print('loaded hdulist... in (s):',f-s)


        for cn in camnums:
            cs='{:02d}'.format(cn)
            cam_dir=os.path.join(output_dir, 'hires_images_cam'+cs)
            print(cam_dir)
        
            if not os.path.lexists(cam_dir):
                os.mkdir(cam_dir)

            export_hdst_filters(hdu_obj,cn,sim_label+'_'+time_label,outdir=cam_dir,nonscatter=nonscatter)
        

            
    
    return
