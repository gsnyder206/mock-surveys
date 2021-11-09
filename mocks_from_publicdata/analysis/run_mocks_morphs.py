import numpy as np
import numpy.random as random
import astropy
import astropy.io.fits as fits
from astropy.convolution import Gaussian2DKernel
import photutils
import statmorph
from statmorph.utils.image_diagnostics import make_figure
import matplotlib
import matplotlib.pyplot as pyplot
from PIL import Image
from astropy.convolution import *
from astropy.stats import gaussian_fwhm_to_sigma
import scipy.ndimage
import scipy as sp
import importlib as imp
import photutils
from photutils import aperture_photometry
from photutils import CircularAperture
from photutils.segmentation import SourceCatalog
from photutils.segmentation import deblend_sources
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy import units as u
import astropy.nddata
from astropy.nddata import Cutout2D
import os
import pickle


def do_photutils_stuff(image,nlevels=32,contrast=0.001,npixels=10):


    threshold = photutils.detect_threshold(image, 1.5)
    kernel_sigma = 3.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))  # FWHM = 3
    kernel = Gaussian2DKernel(kernel_sigma, x_size=9, y_size=9)
    kernel.normalize()
    segm = photutils.detect_sources(image, threshold, npixels, kernel=kernel)

    cat = SourceCatalog(image, segm)
    tbl = cat.to_table()

    #need a deblending step here!
    segm_deblend = deblend_sources(image, segm, kernel=kernel, npixels=npixels,nlevels=nlevels, contrast=contrast)

    db_cat = SourceCatalog(image,segm_deblend)
    db_tbl=db_cat.to_table()

    return segm_deblend,db_tbl,db_cat

def add_simple_noise_extractedsn(in_hdu,radius_arcsec=0.5,extractedsn=300):


    image_in=in_hdu.data
    header_in=in_hdu.header

    #this is the approximate equation for sb magnitude input limit
    #sigma_njy=(2.0**(-0.5))*((1.0e9)*(3631.0/5.0)*10.0**(-0.4*sb_maglim))*header_in['PIXSIZE']*(3.0*header_in['FWHM'])

    #get flux in aperture
    #can we guarantee that the galaxy is in the middle?
    npix=image_in.shape[0]
    ci=np.float32(npix)/2

    radius_pixels = radius_arcsec/in_hdu.header['PIXSCALE']

    positions = [(ci, ci),]
    aperture = CircularAperture(positions, r=radius_pixels)
    phot_table = aperture_photometry(image_in, aperture)
    flux_aperture=phot_table['aperture_sum'][0]

    #get npix in aperture
    area_pixels=np.pi*radius_pixels**2



    #convert to pixel noise level
    sigma_njy=flux_aperture/(extractedsn*(area_pixels)**0.5)



    noise_image = sigma_njy*np.random.randn(npix,npix)

    image_out=image_in + noise_image

    hdu_out = fits.ImageHDU(image_out,header=header_in)
    hdu_out.header['EXTNAME']='MockImage_SN'
    hdu_out.header['EXTSN']=(extractedsn,'rough extracted S/N ratio')
    hdu_out.header['APERFLUX']=flux_aperture
    hdu_out.header['APERRAD']=radius_pixels
    hdu_out.header['RMSNOISE']=(sigma_njy,'nanojanskies')

    return hdu_out

def convolve_with_fwhm_and_rebin(in_hdu, fwhm_arcsec=0.10, desired_pix_side_arcsec=None):

    #load image data and metadata
    image_in=in_hdu.data
    header_in=in_hdu.header

    extlabel=header_in['EXTNAME'].split('_')[-1]

    pixel_size_arcsec=header_in['pixscale']

    sigma_arcsec=fwhm_arcsec*gaussian_fwhm_to_sigma
    sigma_pixels=sigma_arcsec/pixel_size_arcsec

    image_out=sp.ndimage.filters.gaussian_filter(image_in,sigma_pixels,mode='nearest')


    if desired_pix_side_arcsec is not None:
        np_orig=image_out.shape[0]

        np_new=np_orig*pixel_size_arcsec/desired_pix_side_arcsec
        np_new_int=np.int64(np_new)

        orig_fov = np_orig*pixel_size_arcsec
        new_fov = np_new_int*desired_pix_side_arcsec
        diff = (orig_fov-new_fov)/2.0

        box_arcsec=(diff,diff,orig_fov-diff,orig_fov-diff)
        box=(diff/pixel_size_arcsec,diff/pixel_size_arcsec,(orig_fov-diff)/pixel_size_arcsec,(orig_fov-diff)/pixel_size_arcsec)

        #multiply by pixel scale ratio squared in order to preserve total flux
        rebinned_image=np.asarray(Image.fromarray(image_out).resize(size=(np_new_int, np_new_int),box=box))*(desired_pix_side_arcsec/pixel_size_arcsec)**2

        final_image=rebinned_image
        out_pix_arcsec = desired_pix_side_arcsec
    else:
        final_image=image_out
        out_pix_arcsec = pixel_size_arcsec

    hdu_out = fits.ImageHDU(final_image,header=header_in)
    hdu_out.header['FWHM']=(fwhm_arcsec,'arcsec')
    hdu_out.header['SIGMA']=(sigma_arcsec,'arcsec')
    hdu_out.header['EXTNAME']='MockData_PSF_'+extlabel
    hdu_out.header['PIXSCALE']=(out_pix_arcsec,'arcsec')
    hdu_out.header['PIXORIG']=(pixel_size_arcsec,'arcsec, pristine image')
    hdu_out.header['IN_EXT']=header_in['EXTNAME']

    return hdu_out

def source_match(test_source,mockcatalog,image,segm,image_wcs,zoom=3,show=True):

    source_sc = image_wcs.pixel_to_world(test_source.xcentroid,test_source.ycentroid)
    dx=test_source.bbox.ixmax - test_source.bbox.ixmin
    dy=test_source.bbox.iymax - test_source.bbox.iymin

    cut=Cutout2D(image,source_sc,(dy*zoom,dx*zoom),wcs=image_wcs,mode='strict')


    catpix_x,catpix_y=cut.wcs.world_to_pixel(SkyCoord(mockcatalog.data['RA degree']*u.deg,u.deg*mockcatalog.data['DEC degree']))

    catpix_x_im,catpix_y_im=image_wcs.world_to_pixel(SkyCoord(mockcatalog.data['RA degree']*u.deg,u.deg*mockcatalog.data['DEC degree']))


    psep = ((catpix_x_im-test_source.xcentroid)**2 + (catpix_y_im-test_source.ycentroid)**2)**0.5

    fratio= np.abs(test_source.segment_flux - mockcatalog.data['total_quant'])/test_source.segment_flux

    candidates = np.logical_and(psep < 25, fratio < 0.5)

    #if 1, assume decent match and proceed.
    if np.sum(candidates)==1:
        catalog_match=mockcatalog.data[candidates]
        colnames=mockcatalog.data.names
        catdict={}
        catdict['cat_i']=np.where(candidates==True)[0]

        for name,datapoint in zip(colnames,catalog_match[0]):
            catdict[name]=datapoint


        catdict['x_pix']=catpix_x_im[candidates][0]
        catdict['y_pix']=catpix_y_im[candidates][0]

        seg_cut = Cutout2D(segm.data, (test_source.xcentroid,test_source.ycentroid), (dy*zoom,dx*zoom), mode='strict')

        catdict['seg_cutout']=seg_cut.data

        return catdict
    else:
        return None

    #cpx=catpix_x - cut.
    #cpy=catpix_y - cut.

    if show is True:
        f=pyplot.figure(figsize=(8,4),dpi=600)
        pyplot.subplots_adjust(wspace=0.0,hspace=0.0,top=1.0,right=1.0,left=0.00,bottom=0.00)

        ax=f.add_subplot(141)
        ax.set_xticks([]) ; ax.set_yticks([])

        ax.imshow(test_source.data,interpolation='nearest',aspect='auto',origin='lower')

        ax=f.add_subplot(142)
        ax.set_xticks([]) ; ax.set_yticks([])

        ax.imshow(test_source.segment,interpolation='nearest',aspect='auto',origin='lower')

        ax=f.add_subplot(143)
        ax.set_xticks([]) ; ax.set_yticks([])


        ax.imshow(cut.data,interpolation='nearest',aspect='auto')

        ax.plot(catpix_x[candidates],catpix_y[candidates],'+',color='red',markersize=10)

        ax.set_xlim(0,dx*zoom) ; ax.set_ylim(0,dy*zoom)

        ax.set_xticks([]) ; ax.set_yticks([])



        ax=f.add_subplot(144)
        ax.set_xticks([]) ; ax.set_yticks([])



    return candidates

def get_catline_header():
    names=('snapnum','subfindid','redshift','total_quant_match','mstar','mgas','msub','mbh','mbar','sfr','bhar','x_pix','y_pix','total_quant',
        'photutils_label','photutils_xcen','photutils_ycen','segment_flux',
        'sm_flag','sm_flag_sersic','sm_flag_catastrophic',
        'sm_gini','sm_m20','sm_bulge','sm_merger','sm_sn_per_pixel','sm_c','sm_a','sm_rhalf_e','sm_rpetro_e','sm_flux_e',
        'sm_xc_cen','sm_yc_cen','sm_sersic_n','sm_sersic_rhalf','sm_shape_asym','sm_m','sm_i','sm_d','sm_label')

    no=('snapnum','subfindid','redshift','total_quant','mstar','mgas','msub','mbh','mbar','sfr','bhar','x_pix','y_pix',
        'photutils_label','photutils_xcen','photutils_ycen',
        'sm_flag','sm_flag_sersic',
        'sm_gini','sm_m20','sm_bulge','sm_merger','sm_sn_per_pixel','sm_c','sm_a','sm_rhalf_e','sm_rpetro_e','sm_flux_e',
        'sm_xc_cen','sm_yc_cen','sm_sersic_n','sm_sersic_rhalf','sm_shape_asym','sm_m','sm_i','sm_d','sm_label')
    header='{} {} {} {} {} {} {} {} {} {} {} {} {} {} '\
            '{} {} {} {} '\
            '{} {} {} {} {} {} {} '\
            '{} {} {} {} {} {} {} '\
            '{} {} {} {} {} {} {} {}\n'.format(*names)
    return header,names,no

def get_catline(so,cdict,sm,this_cat,filter_pu_flux):
    this_cat_dict={}
    this_cat_row=this_cat.data[cdict['cat_i']]
    for dp,cn in zip(this_cat_row[0],this_cat.data.names):
        this_cat_dict[cn]=dp

    #print(this_cat_dict)
    stuff_tuple=(cdict['Snapshot number'],
                cdict['Subhalo index'],
                cdict['True z'],
                cdict['total_quant'],
                cdict['Stellar mass w2sr'],
                cdict['Total gas mass w2sr'],
                cdict['Total subhalo mass'],
                cdict['Total BH mass w2sr'],
                cdict['Total baryon mass w2sr'],
                cdict['SFR w2sr'],
                cdict['Total BH accretion rate'],
                cdict['x_pix'],
                cdict['y_pix'],
                this_cat_dict['total_quant'],
                so.label,
                so.xcentroid,
                so.ycentroid,
                filter_pu_flux,
                sm.flag,
                sm.flag_sersic,
                sm.flag_catastrophic,
                sm.gini,
                sm.m20,
                sm.gini_m20_bulge,
                sm.gini_m20_merger,
                sm.sn_per_pixel,
                sm.concentration,
                sm.asymmetry,
                sm.rhalf_ellip,
                sm.rpetro_ellip,
                sm.flux_ellip,
                sm.xc_centroid,
                sm.yc_centroid,
                sm.sersic_n,
                sm.sersic_rhalf,
                sm.shape_asymmetry,
                sm.multimode,
                sm.intensity,
                sm.deviation,
                sm.label)

    line='{:10d} {:10d} {:15.8f} '\
        '{:12.8f} {:12.8} {:12.8} {:12.8} {:12.8} {:12.8} {:12.8f} {:12.8f} '\
        '{:10.6f} {:10.6f} {:12.8f} '\
        '{:10d} {:10.6f} {:10.6f} {:12.8f} '\
        '{:10d} {:10d} {:10d} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} '\
        '{:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} '\
        '{:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10d}\n'.format(*stuff_tuple)

    return line

def make_sm_fig(source_morph,thisdir,filelabel):
    fig = make_figure(source_morph)
    fig.savefig(os.path.join(thisdir,'statmorph_'+filelabel+'.png'),dpi=150)
    pyplot.close(fig)
    return


def run_morphs_calculations(im435,im814,im200,im277,im356,im_mstar,
                            cat_435,cat_814,cat_200,cat_277,cat_356,cat_mstar,
                            wcs_200,segm_1,tbl_1,scat_1,label="test",min_mag=24.0,max_mag=29.0,limit=10,starti=0,
                            outdir_root='/astro/snyder_lab/MockSurveys/MockLightcones/morph_analysis/'):
    im_200_1=1.0*im200

    #segm_1,tbl_1,scat_1 = do_photutils_stuff(im_200_1)



    totals=tbl_1['segment_flux']
    labels = tbl_1['label']
    mags = -2.5*np.log10(1.0e-9*totals/3631.0)
    si=np.logical_and(mags < max_mag, mags > min_mag)

    segm_1_all = segm_1.copy()
    segm_1.keep_labels(labels[si])
    print(len(tbl_1),np.sum(si))

    n=1000
    f=pyplot.figure(figsize=(20,10),dpi=200)
    pyplot.subplots_adjust(wspace=0.0,hspace=0.0,top=1.0,right=1.0,left=0.00,bottom=0.00)
    ax=f.add_subplot(121)
    ax.imshow(np.log10(im_200_1[0:n,0:n]))
    ax.set_xticks([]) ; ax.set_yticks([])

    catpix_x,catpix_y=wcs_200.world_to_pixel(SkyCoord(cat_200.data['RA degree']*u.deg,u.deg*cat_200.data['DEC degree']))
    ax.plot(catpix_x,catpix_y,'+',color='red')
    ax.set_xlim(0,n) ; ax.set_ylim(0,n)

    ax=f.add_subplot(122)
    ax.set_xticks([]) ; ax.set_yticks([])
    cmap = segm_1.make_cmap(seed=123)
    pyplot.imshow(segm_1.data[0:n,0:n],cmap=cmap)
    ax.plot(catpix_x,catpix_y,'+',color='red')
    ax.set_xlim(0,n) ; ax.set_ylim(0,n)
    f.savefig('image_segments_example_'+label+'.png',dpi=200)
    pyplot.close(f)




    outputs_dir=os.path.join(outdir_root,'morphology_diagnostics_'+label)
    f200_dir=os.path.join(outputs_dir,'f200w')
    f435_dir=os.path.join(outputs_dir,'f435w')
    f814_dir=os.path.join(outputs_dir,'f814w')
    f277_dir=os.path.join(outputs_dir,'f277w')
    f356_dir=os.path.join(outputs_dir,'f356w')
    mstar_dir=os.path.join(outputs_dir,'mstar')

    if not os.path.lexists(outputs_dir):
        os.makedirs(outputs_dir)
    if not os.path.lexists(f200_dir):
        os.makedirs(f200_dir)
    if not os.path.lexists(f435_dir):
        os.makedirs(f435_dir)
    if not os.path.lexists(f814_dir):
        os.makedirs(f814_dir)
    if not os.path.lexists(f277_dir):
        os.makedirs(f277_dir)
    if not os.path.lexists(f356_dir):
        os.makedirs(f356_dir)
    if not os.path.lexists(mstar_dir):
        os.makedirs(mstar_dir)



    morphcat_200 = open(os.path.join(f200_dir,'morphcat_'+label+'_f200w.txt'),'w')
    morphcat_435 = open(os.path.join(f435_dir,'morphcat_'+label+'_f435w.txt'),'w')
    morphcat_814 = open(os.path.join(f814_dir,'morphcat_'+label+'_f814w.txt'),'w')
    morphcat_277 = open(os.path.join(f277_dir,'morphcat_'+label+'_f277w.txt'),'w')
    morphcat_356 = open(os.path.join(f356_dir,'morphcat_'+label+'_f356w.txt'),'w')
    morphcat_mstar=open(os.path.join(mstar_dir,'morphcat_'+label+'_mstar.txt'),'w')

    header,names,no=get_catline_header()
    morphcat_200.write(header)
    morphcat_435.write(header)
    morphcat_814.write(header)
    morphcat_277.write(header)
    morphcat_356.write(header)
    morphcat_mstar.write(header)



    if limit is not None:
        iterable_thing = scat_1[si][starti:starti+limit]
    else:
        iterable_thing = scat_1[si][starti:]

    N_segments = len(iterable_thing)
    zoom=3
    for i,test_source in enumerate(iterable_thing):
        cdict=source_match(test_source,cat_200,im_200_1,segm_1,wcs_200,show=False)
        if i % 100 == 0:
            print(i, ' segments completed of ', N_segments)


        if cdict!=None:

            #may have to replace this with using cutouts?
            #seg_copy=segm_1.copy()
            #seg_copy.keep_label(test_source.label)
            seg_label = cdict['seg_cutout']==test_source.label
            dx=test_source.bbox.ixmax - test_source.bbox.ixmin
            dy=test_source.bbox.iymax - test_source.bbox.iymin
            #catdict['im_cutout']=cut.data
            f200w_cut = Cutout2D(im_200_1, (test_source.xcentroid,test_source.ycentroid), (dy*zoom,dx*zoom), mode='strict')
            f435w_cut = Cutout2D(im_435_1, (test_source.xcentroid,test_source.ycentroid), (dy*zoom,dx*zoom), mode='strict')
            f814w_cut = Cutout2D(im_814_1, (test_source.xcentroid,test_source.ycentroid), (dy*zoom,dx*zoom), mode='strict')
            f356w_cut = Cutout2D(im_356_1, (test_source.xcentroid,test_source.ycentroid), (dy*zoom,dx*zoom), mode='strict')
            f277w_cut = Cutout2D(im_277_1, (test_source.xcentroid,test_source.ycentroid), (dy*zoom,dx*zoom), mode='strict')
            mstar_cut = Cutout2D(im_mstar_1, (test_source.xcentroid,test_source.ycentroid), (dy*zoom,dx*zoom), mode='strict')

            use_seg_label = np.logical_or(cdict['seg_cutout'] == 0 , cdict['seg_cutout']==test_source.label)

            #have to mask out back/fore ground sources.
            source_morph_200 = statmorph.source_morphology(np.where(use_seg_label,f200w_cut.data,0.0), np.where(seg_label,cdict['seg_cutout'],0), gain=5)
            source_morph_435 = statmorph.source_morphology(np.where(use_seg_label,f435w_cut.data,0.0), np.where(seg_label,cdict['seg_cutout'],0), gain=5)
            source_morph_814 = statmorph.source_morphology(np.where(use_seg_label,f814w_cut.data,0.0), np.where(seg_label,cdict['seg_cutout'],0), gain=5)
            source_morph_356 = statmorph.source_morphology(np.where(use_seg_label,f356w_cut.data,0.0), np.where(seg_label,cdict['seg_cutout'],0), gain=5)
            source_morph_277 = statmorph.source_morphology(np.where(use_seg_label,f277w_cut.data,0.0), np.where(seg_label,cdict['seg_cutout'],0), gain=5)
            source_morph_mstar = statmorph.source_morphology(np.where(use_seg_label,mstar_cut.data,0.0), np.where(seg_label,cdict['seg_cutout'],0), gain=5)

            seg_flux_200=np.sum(np.where(seg_label,f200w_cut.data,0.0))
            seg_flux_435=np.sum(np.where(seg_label,f435w_cut.data,0.0))
            seg_flux_814=np.sum(np.where(seg_label,f814w_cut.data,0.0))
            seg_flux_277=np.sum(np.where(seg_label,f277w_cut.data,0.0))
            seg_flux_356=np.sum(np.where(seg_label,f356w_cut.data,0.0))
            seg_flux_mstar=np.sum(np.where(seg_label,mstar_cut.data,0.0))


            #source_morph_435 = statmorph.source_morphology(np.where(segm_1_all.data==seg_copy.data,im_435_1,0.0), seg_copy.data, gain=5)
            #source_morph_814 = statmorph.source_morphology(np.where(segm_1_all.data==seg_copy.data,im_814_1,0.0), seg_copy.data, gain=5)
            #source_morph_356 = statmorph.source_morphology(np.where(segm_1_all.data==seg_copy.data,im_356_1,0.0), seg_copy.data, gain=5)
            #source_morph_277 = statmorph.source_morphology(np.where(segm_1_all.data==seg_copy.data,im_277_1,0.0), seg_copy.data, gain=5)
            #source_morph_mstar = statmorph.source_morphology(np.where(segm_1_all.data==seg_copy.data,im_mstar_1,0.0), seg_copy.data, gain=5)


            filelabel=str(test_source.label)+'_'+str(cdict['Snapshot number'])+'_'+str(cdict['Subhalo index'])

            try:
                make_sm_fig(source_morph_200[0],f200_dir,filelabel)
            except AssertionError as AE:
                pass
            try:
                make_sm_fig(source_morph_435[0],f435_dir,filelabel)
            except AssertionError as AE:
                pass
            try:
                make_sm_fig(source_morph_814[0],f814_dir,filelabel)
            except AssertionError as AE:
                pass
            try:
                make_sm_fig(source_morph_277[0],f277_dir,filelabel)
            except AssertionError as AE:
                pass
            try:
                make_sm_fig(source_morph_356[0],f356_dir,filelabel)
            except AssertionError as AE:
                pass
            try:
                make_sm_fig(source_morph_mstar[0],mstar_dir,filelabel)
            except AssertionError as AE:
                pass

            catline_200=get_catline(test_source,cdict,source_morph_200[0],cat_200,seg_flux_200)
            morphcat_200.write(catline_200)
            catline_435=get_catline(test_source,cdict,source_morph_435[0],cat_435,seg_flux_435)
            morphcat_435.write(catline_435)
            catline_814=get_catline(test_source,cdict,source_morph_814[0],cat_814,seg_flux_814)
            morphcat_814.write(catline_814)
            catline_277=get_catline(test_source,cdict,source_morph_277[0],cat_277,seg_flux_277)
            morphcat_277.write(catline_277)
            catline_356=get_catline(test_source,cdict,source_morph_356[0],cat_356,seg_flux_356)
            morphcat_356.write(catline_356)
            catline_mstar=get_catline(test_source,cdict,source_morph_mstar[0],cat_mstar,seg_flux_mstar)
            morphcat_mstar.write(catline_mstar)

            morphcat_200.flush()
            morphcat_435.flush()
            morphcat_814.flush()
            morphcat_277.flush()
            morphcat_356.flush()
            morphcat_mstar.flush()

    morphcat_200.close()
    morphcat_435.close()
    morphcat_814.close()
    morphcat_277.close()
    morphcat_356.close()
    morphcat_mstar.close()
    return

if __name__=="__main__":
    #make some intermediate things

    #im_base='../outputs/TNG50-1_12_11_xyz_'
    im_base='../outputs/TNG100-1_7_6_xyz/TNG100-1_xyz_'

    #label='tng50_test'
    #label='tng100_test'
    label="tng100_run_2021_11_03"

    limit=None
    #limit=20

    cat_435=fits.open(im_base+'F435W.fits')[2]
    cat_814=fits.open(im_base+'F814W.fits')[2]
    cat_200=fits.open(im_base+'F200W.fits')[2]
    cat_277=fits.open(im_base+'F277W.fits')[2]
    cat_356=fits.open(im_base+'F356W.fits')[2]
    cat_mstar=fits.open(im_base+'mstar.fits')[2]

    s1 = 0.02 ; s2 = 0.20 ; s3 = 2.0 #noise levels in njy


    tmp_dir=os.path.join(os.path.dirname(im_base),'tmp_mocks')
    if not os.path.lexists(tmp_dir):
        os.makedirs(tmp_dir)

    try:
        im_435_1=fits.open(os.path.join(tmp_dir,'tmp_f435w_1.fits'))[0].data
        im_814_1=fits.open(os.path.join(tmp_dir,'tmp_f814w_1.fits'))[0].data
        im_200_1=fits.open(os.path.join(tmp_dir,'tmp_f200w_1.fits'))[0].data
        im_277_1=fits.open(os.path.join(tmp_dir,'tmp_f277w_1.fits'))[0].data
        im_356_1=fits.open(os.path.join(tmp_dir,'tmp_f356w_1.fits'))[0].data
        im_mstar_1=fits.open(os.path.join(tmp_dir,'tmp_mstar_1.fits'))[0].data

        ph_200=fits.open(im_base+'F200W.fits')[0]

        #load seg map object stuff
        with open(os.path.join(tmp_dir,'tmp_photutils_results_1.pkl'), 'rb') as inp:
            list_of_seg_stuff = pickle.load(inp)
            segm_1=list_of_seg_stuff[0]
            tbl_1=list_of_seg_stuff[1]
            scat_1=list_of_seg_stuff[2]

    except:
        print('tmp files do not exist, remaking them now... (might be slow)')
        ph_435=fits.open(im_base+'F435W.fits')[0]
        ph_814=fits.open(im_base+'F814W.fits')[0]
        ph_200=fits.open(im_base+'F200W.fits')[0]
        ph_277=fits.open(im_base+'F277W.fits')[0]
        ph_356=fits.open(im_base+'F356W.fits')[0]
        ph_mstar=fits.open(im_base+'mstar.fits')[0]

        psfh_435=convolve_with_fwhm_and_rebin(ph_435,fwhm_arcsec=0.04)
        psfh_814=convolve_with_fwhm_and_rebin(ph_814,fwhm_arcsec=0.08)
        psfh_200=convolve_with_fwhm_and_rebin(ph_200,fwhm_arcsec=0.08)
        psfh_356_03=convolve_with_fwhm_and_rebin(ph_356,fwhm_arcsec=0.12)
        psfh_277_03=convolve_with_fwhm_and_rebin(ph_277,fwhm_arcsec=0.09)
        #psfh_356_06=convolve_with_fwhm_and_rebin(ph_356,fwhm_arcsec=0.12,desired_pix_side_arcsec=0.06)

        psfh_mstar_03=convolve_with_fwhm_and_rebin(ph_mstar,fwhm_arcsec=0.08)

        npix=psfh_200.data.shape[0]

        im_200_1=psfh_200.data+s1*np.random.randn(npix,npix)
        #im_200_2=psfh_200.data+s2*np.random.randn(npix,npix)
        #im_200_3=psfh_200.data+s3*np.random.randn(npix,npix)

        im_435_1=psfh_435.data+s1*np.random.randn(npix,npix)
        im_814_1=psfh_814.data+s1*np.random.randn(npix,npix)
        im_277_1=psfh_277_03.data+s1*np.random.randn(npix,npix)
        im_356_1=psfh_356_03.data+s1*np.random.randn(npix,npix)
        im_mstar_1=psfh_mstar_03.data+0.1*np.random.randn(npix,npix) #made up mstar noise

        fits.PrimaryHDU(im_435_1).writeto(os.path.join(tmp_dir,'tmp_f435w_1.fits'),overwrite=True)
        fits.PrimaryHDU(im_814_1).writeto(os.path.join(tmp_dir,'tmp_f814w_1.fits'),overwrite=True)
        fits.PrimaryHDU(im_200_1).writeto(os.path.join(tmp_dir,'tmp_f200w_1.fits'),overwrite=True)
        fits.PrimaryHDU(im_277_1).writeto(os.path.join(tmp_dir,'tmp_f277w_1.fits'),overwrite=True)
        fits.PrimaryHDU(im_356_1).writeto(os.path.join(tmp_dir,'tmp_f356w_1.fits'),overwrite=True)
        fits.PrimaryHDU(im_mstar_1).writeto(os.path.join(tmp_dir,'tmp_mstar_1.fits'),overwrite=True)

        segm_1,tbl_1,scat_1 = do_photutils_stuff(np.float32(im_200_1))
        with open(os.path.join(tmp_dir,'tmp_photutils_results_1.pkl'), 'wb') as outp:
            pickle.dump([segm_1,tbl_1,scat_1], outp)



    wcs_200 = wcs.WCS(ph_200)


    run_morphs_calculations(im_435_1,im_814_1,im_200_1,im_277_1,im_356_1,im_mstar_1,
                            cat_435,cat_814,cat_200,cat_277,cat_356,cat_mstar,
                            wcs_200,segm_1,tbl_1,scat_1,label=label,min_mag=20.0,max_mag=30.0,limit=limit,starti=47000)
