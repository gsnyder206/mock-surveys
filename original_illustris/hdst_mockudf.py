import math
import string
import sys
import struct
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.colors as pycolors
import matplotlib.cm as cm
import matplotlib.patches as patches
import numpy as np
import cPickle
import asciitable
import scipy.ndimage
import scipy.stats as ss
import scipy.signal
import scipy as sp
import scipy.odr as odr
import astropy.io.fits as pyfits
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

sq_arcsec_per_sr = 42545170296.0
c = 3.0e8

def render_only(outfile='HDUDF_v1.pdf',hst_only=False,maglim=28,label='_SB28_',unit='nJySr'):

    print "reading b"
    b=pyfits.open('hudf_F606W_Jy.fits')[0].data #hdudf_v.final_array
    print "reading g"
    g=pyfits.open('hudf_F850LP_Jy.fits')[0].data #hdudf_z.final_array
    print "reading r"
    r=pyfits.open('hudf_F160W_Jy.fits')[0].data #hdudf_h.final_array
    pixel_arcsec = pyfits.open('hudf_F160W_Jy.fits')[0].header['PIXSCALE']
    if unit=='Jy':
        conv = (1.0e9)*(1.0/pixel_arcsec**2)*sq_arcsec_per_sr

    #res = render_hdudf(b*conv,g*conv,r*conv,'HUDF'+label+'_v1.pdf',pixel_arcsec=pixel_arcsec,FWHM_arcsec_b=0.10,FWHM_arcsec_g=0.15,FWHM_arcsec_r=0.20,convolve=True,dpi=600,maglim=maglim)
    #res = render_hdudf(b,g,r,'HUDF'+label+'small_v1.jpg',pixel_arcsec=pixel_arcsec,FWHM_arcsec_b=0.10,FWHM_arcsec_g=0.15,FWHM_arcsec_r=0.20,convolve=True,dpi=60,maglim=maglim)
    res = render_hdudf(b*conv,g*conv,r*conv,'HUDF'+label+'big_v4.jpg',pixel_arcsec=pixel_arcsec,FWHM_arcsec_b=0.10,FWHM_arcsec_g=0.15,FWHM_arcsec_r=0.20,convolve=True,dpi=1200,maglim=maglim)
    res = render_hdudf(b*conv,g*conv,r*conv,'HUDF'+label+'small_v4.jpg',pixel_arcsec=pixel_arcsec,FWHM_arcsec_b=0.10,FWHM_arcsec_g=0.15,FWHM_arcsec_r=0.20,convolve=True,dpi=60,maglim=maglim)
    if hst_only==True:
        return

    print "reading b"
    b=pyfits.open('hdudf_6mas_F606W_Jy.fits')[0].data #hdudf_v.final_array
    print "reading g"
    g=pyfits.open('hdudf_6mas_F850LP_Jy.fits')[0].data#hdudf_z.final_array
    print "reading r"
    r=pyfits.open('hdudf_6mas_F160W_Jy.fits')[0].data#hdudf_h.final_array
    pixel_arcsec = pyfits.open('hdudf_6mas_F160W_Jy.fits')[0].header['PIXSCALE']
    if unit=='Jy':
        conv = (1.0e9)*(1.0/pixel_arcsec**2)*sq_arcsec_per_sr

    #assume trying for 8m telescope
    #res = render_hdudf(b*conv,g*conv,r*conv,'HDUDF'+label+'_v1.pdf',pixel_arcsec=pixel_arcsec,FWHM_arcsec_b=0.017,FWHM_arcsec_g=0.025,FWHM_arcsec_r=0.050,convolve=True,dpi=2000,maglim=maglim)

    #settings for 12m
    res = render_hdudf(b*conv,g*conv,r*conv,'HDUDF'+label+'big_v4.jpg',pixel_arcsec=pixel_arcsec,FWHM_arcsec_b=0.012,FWHM_arcsec_g=0.018,FWHM_arcsec_r=0.032,convolve=True,dpi=1200,maglim=maglim)


    return



def render_hdudf(b,g,r,filename,pixel_arcsec=0.006,FWHM_arcsec_b=0.012,FWHM_arcsec_g=0.015,FWHM_arcsec_r=0.025,convolve=True,dpi=2000,maglim=28.0):
    #maglim in mag/arcsec^2

    redfact = 1.5*(0.60/1.60)**(1)
    greenfact = 0.9*(0.60/0.85)**(1)
    bluefact = 1.2
    efflams = [1.60,1.25,0.90,0.775,0.606,0.435,0.814,1.05,1.40]

    alph=7.0
    Q = 5.0

    target_ratio = 10.0**(-0.4*(27.0-maglim))
    fluxscale = target_ratio*1.0e-14
    pixel_Sr = (pixel_arcsec**2)/sq_arcsec_per_sr


    #new version, already in nJy/Sr
    to_nJy_per_Sr_b = 1#(1.0e9)*(1.0e14)*(efflams[4]**2)/c   #((pixscale/206265.0)^2)*
    to_nJy_per_Sr_g = 1#(1.0e9)*(1.0e14)*(efflams[2]**2)/c
    to_nJy_per_Sr_r = 1#(1.0e9)*(1.0e14)*(efflams[0]**2)/c

    #b_nJySr = to_nJy_per_Sr_b*b
    #g_nJySr = to_nJy_per_Sr_g*g
    #r_nJySr = to_nJy_per_Sr_r*r


    sigma_pixels_b = FWHM_arcsec_b/pixel_arcsec/2.355
    sigma_pixels_g = FWHM_arcsec_g/pixel_arcsec/2.355
    sigma_pixels_r = FWHM_arcsec_r/pixel_arcsec/2.355

    print "sigma pixels: ", sigma_pixels_g
    if convolve==True:
        print "convolving images"
        b = scipy.ndimage.filters.gaussian_filter(b,sigma_pixels_b)
        g = scipy.ndimage.filters.gaussian_filter(g,sigma_pixels_g)
        r = scipy.ndimage.filters.gaussian_filter(r,sigma_pixels_r)



    sigma_nJy = 0.3*(2.0**(-0.5))*((1.0e9)*(3631.0/5.0)*10.0**(-0.4*maglim))*pixel_arcsec*(3.0*FWHM_arcsec_g)

    print "adding noise, in nJy/Sr: ", sigma_nJy/pixel_Sr
    Npix = b.shape[0]
    b = (b*to_nJy_per_Sr_b + np.random.randn(Npix,Npix)*sigma_nJy/pixel_Sr)
    g = (g*to_nJy_per_Sr_g + np.random.randn(Npix,Npix)*sigma_nJy/pixel_Sr)
    r = (r*to_nJy_per_Sr_r + np.random.randn(Npix,Npix)*sigma_nJy/pixel_Sr)



    print "preparing color image"
    rgbdata = make_color_image.make_interactive_light_nasa(b*fluxscale*bluefact,g*fluxscale*greenfact,r*fluxscale*redfact,alph,Q)

    print rgbdata.shape
    print "preparing figure"
    f9 = pyplot.figure(figsize=(12.0,12.0), dpi=dpi)
    pyplot.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0,wspace=0.0,hspace=0.0)

    axi = pyplot.axes([0.0,0.0,1.0,1.0],frameon=True,axisbg='black')
    axi.set_xticks([]) ; axi.set_yticks([])
    print "rendering color image"
    axi.imshow(rgbdata,interpolation='nearest',origin='upper',extent=[-1,1,-1,1])  

    print "saving color image"
    #pyplot.rcParams['pdf.compression'] = 1
    f9.savefig(filename,dpi=dpi,quality=90,bbox_inches='tight',pad_inches=0.0)
    pyplot.close(f9)
    #pyplot.rcdefaults()



    return

class mock_hdudf:
    def __init__(self,Npix,Pix_arcsec,blank_array,filter_string,simdata,Narcmin,eff_lambda_microns,maglim,fwhm,req_filters=[]):
        self.Npix=Npix
        self.Pix_arcsec=Pix_arcsec
        self.fov_arcmin = Narcmin
        self.blank_array=blank_array*1.0
        self.final_array=blank_array*1.0
        self.filter_string=filter_string
        self.simdata=simdata
        self.image_files=[]
        self.x_array=[]
        self.y_array=[]
        self.N_inf=[]
        self.eff_lambda_microns = eff_lambda_microns
        self.maglim = 28.0
        self.FWHM_arcsec = fwhm
        self.req_filters=req_filters
        self.mstar_list = []
        self.redshift_list = []

    def find_image(self,mstar,redshift,sfr,seed,xpix,ypix,hmag):
        sim_simname = self.simdata['col1']
        sim_expfact = self.simdata['col2']
        sim_sfr = self.simdata['col54']
        sim_mstar = self.simdata['col56']
        sim_redshift = 1.0/sim_expfact - 1.0
        metalmass = self.simdata['col53']
        sim_res_pc = self.simdata['col62']
        sim_string = self.simdata['col60']

        simage_loc = '/Users/gsnyder/Documents/Projects/HydroART_Morphology/Hyades_Data/images_rsync/'

        self.mstar_list.append(mstar)
        self.redshift_list.append(redshift)

        adjust_size=False

        print " "
        print "Searching for simulation with mstar,z,seed : ", mstar, redshift, seed
        wide_i = np.where(np.logical_and(np.logical_and(np.abs(sim_redshift-redshift)<0.3,np.abs(np.log10(sim_mstar)-mstar)<0.1),sim_sfr > -1))[0]
        Nwi = wide_i.shape[0]
        if Nwi==0:
            wide_i = np.where(np.logical_and(np.logical_and(np.abs(sim_redshift-redshift)<0.5,np.abs(np.log10(sim_mstar)-mstar)<0.4),sim_sfr > -1))[0]
            Nwi = wide_i.shape[0]
        if Nwi==0 and (mstar < 7.1):
            print "  Can't find good sim, adjusting image parameters to get low mass things "
            wide_i = np.where(np.abs(sim_redshift-redshift)<0.3)[0]  #wide_i is a z range
            llmi = np.argmin(np.log10(sim_mstar[wide_i]))  #the lowest mass in this z range
            wlmi = np.where(np.abs(np.log10(sim_mstar[wide_i]) - np.log10(sim_mstar[wide_i[llmi]])) < 0.3)[0] #search within 0.3 dex of lowest available sims
            print "   ", wide_i.shape, llmi, wlmi.shape
            wide_i = wide_i[wlmi]
            Nwi = wide_i.shape[0]
            print "   ", Nwi
            adjust_size=True

        #assert(wide_i.shape[0] > 0)

        if Nwi==0:
            print "    Could not find roughly appropriate simulation for mstar,z: ", mstar, redshift
            print " "
            self.image_files.append('')
            return 0#np.zeros(shape=(600,600)), -1

        print "    Found N candidates: ", wide_i.shape


        np.random.seed(seed)

        #choose random example and camera
        rps = np.random.random_integers(0,Nwi-1,1)[0]
        cn = str(np.random.random_integers(5,8,1)[0])

        prefix = os.path.basename(sim_string[wide_i[rps]])
        sim_realmstar = np.log10(sim_mstar[wide_i[rps]]) #we picked a sim with this log mstar
        mstar_factor = sim_realmstar - mstar  

        rad_factor = 1.0
        lum_factor = 1.0

        if adjust_size==True:
            rad_factor = 10.0**(mstar_factor*0.5) #must **shrink** images by this factor, total flux by mstar factor
            lum_factor = 10.0**(mstar_factor)
    

        print ">>>FACTORS<<<   ", prefix, sim_realmstar, mstar_factor, rad_factor, lum_factor

        im_folder = simage_loc + prefix +'_skipir/images'

        im_file = os.path.join(im_folder, prefix+'_skipir_CAMERA'+cn+'-BROADBAND_'+self.filter_string+'_simulation.fits')
        cn_file = os.path.join(im_folder, prefix+'_skipir_CAMERA'+cn+'-BROADBAND_'+self.filter_string+'_candelized_noise.fits')
        req1 = os.path.join(im_folder, prefix+'_skipir_CAMERA'+cn+'-BROADBAND_'+self.req_filters[0]+'_simulation.fits')
        req2 = os.path.join(im_folder, prefix+'_skipir_CAMERA'+cn+'-BROADBAND_'+self.req_filters[1]+'_simulation.fits')
        req3 = os.path.join(im_folder, prefix+'_skipir_CAMERA'+cn+'-BROADBAND_'+self.req_filters[2]+'_simulation.fits')
        req4 = os.path.join(im_folder, prefix+'_skipir_CAMERA'+cn+'-BROADBAND_'+self.req_filters[3]+'_simulation.fits')
        req5 = os.path.join(im_folder, prefix+'_skipir_CAMERA'+cn+'-BROADBAND_'+self.req_filters[4]+'_simulation.fits')


        ## Actually, probably want to keep trying some possible galaxies/files...
        is_file = os.path.lexists(im_file) and os.path.lexists(cn_file) and os.path.lexists(req1) and os.path.lexists(req2) and os.path.lexists(req3) and os.path.lexists(req4) and os.path.lexists(req5)
        #is_file = os.path.lexists(im_file) and os.path.lexists(cn_file) #and os.path.lexists(req1) and os.path.lexists(req2) and os.path.lexists(req3)

        if is_file==False:
            print "    Could not find appropriate files: ", im_file, cn_file
            print " "
            self.image_files.append('')
            return 0 #np.zeros(shape=(600,600)), -1


        self.image_files.append(im_file)

        cn_header = pyfits.open(cn_file)[0].header

        im_hdu = pyfits.open(im_file)[0]
        scalesim = cn_header.get('SCALESIM') #pc/pix
        Ps = cosmocalc.cosmocalc(redshift)['PS_kpc'] #kpc/arcsec
        print "    Simulation pixel size at z: ", scalesim
        print "    Plate scale for z: ", Ps
        print "    Desired Kpc/pix at z: ", Ps*self.Pix_arcsec
        sunrise_image = np.float32(im_hdu.data)  #W/m/m^2/Sr

        Sim_Npix = sunrise_image.shape[0]
        New_Npix = int( Sim_Npix*(scalesim/(1000.0*Ps*self.Pix_arcsec))/rad_factor )  #rad_factor reduces number of pixels (total size) desired
        if New_Npix==0:
            New_Npix=1

        print "    New galaxy pixel count: ", New_Npix
        rebinned_image = congrid.congrid(sunrise_image,(New_Npix,New_Npix)) #/lum_factor  #lum_factor shrinks surface brightness by mass factor... but we're shrinking size first, so effective total flux already adjusted by this; may need to ^^ SB instead???  or fix size adjust SB?
        print "    New galaxy image shape: ", rebinned_image.shape
        print "    New galaxy image max: ", np.max(rebinned_image)
        #finite_bool = np.isfinite(rebinned_image)
        #num_infinite = np.where(finite_bool==False)[0].shape[0]
        #print "    Number of INF pixels: ", num_infinite, prefix
        #self.N_inf.append(num_infinite)

        if xpix==-1:
            xpix = int( (self.Npix-1)*np.random.rand()) #np.random.random_integers(0,self.Npix-1,1)[0]
            ypix = int( (self.Npix-1)*np.random.rand()) #np.random.random_integers(0,self.Npix-1,1)[0]

        self.x_array.append(xpix)
        self.y_array.append(ypix)

        x1_choice = np.asarray([int(xpix-float(New_Npix)/2.0),0])
        x1i = np.argmax(x1_choice)
        x1 = x1_choice[x1i]
        diff=0
        if x1==0:
            diff = x1_choice[1]-x1_choice[0]
        x2_choice = np.asarray([x1 + New_Npix - diff,self.Npix])
        x2i = np.argmin(x2_choice)
        x2 = int(x2_choice[x2i])
        x1sim = abs(np.min(x1_choice))
        x2sim = min(New_Npix,self.Npix-x1)

        y1_choice = np.asarray([int(ypix-float(New_Npix)/2.0),0])
        y1i = np.argmax(y1_choice)
        y1 = y1_choice[y1i]
        diff=0
        if y1==0:
            diff = y1_choice[1]-y1_choice[0]
        y2_choice = np.asarray([y1 + New_Npix - diff,self.Npix])
        y2i = np.argmin(y2_choice)
        y2 = int(y2_choice[y2i])
        y1sim = abs(np.min(y1_choice))
        y2sim = min(New_Npix,self.Npix-y1)


        print "    Placing new image at x,y in x1:x2, y1:y2  from xsim,ysim, ", xpix, ypix, x1,x2,y1,y2, x1sim, x2sim, y1sim, y2sim
        #image_slice = np.zeros_like(self.blank_array)
        print "    done creating image slice"
        #bool_slice = np.int32( np.zeros(shape=(self.Npix,self.Npix)))

        image_cutout = rebinned_image[x1sim:x2sim,y1sim:y2sim]
        print "    New image shape: ", image_cutout.shape



        pixel_Sr = (self.Pix_arcsec**2)/sq_arcsec_per_sr  #pixel area in steradians:  Sr/pixel


        to_nJy_per_Sr = (1.0e9)*(1.0e14)*(self.eff_lambda_microns**2)/c   #((pixscale/206265.0)^2)*
        #sigma_nJy = 0.3*(2.0**(-0.5))*((1.0e9)*(3631.0/5.0)*10.0**(-0.4*self.maglim))*self.Pix_arcsec*(3.0*self.FWHM_arcsec)
        to_Jy_per_pix = to_nJy_per_Sr*(1.0e-9)*pixel_Sr
        

        #b = b*(to_nJy_per_Sr_b*fluxscale*bluefact)   # + np.random.randn(Npix,Npix)*sigma_nJy/pixel_Sr
        image_cutout = image_cutout*to_Jy_per_pix #image_cutout*to_nJy_per_Sr


        #image_slice[x1:x2,y1:y2] = image_cutout*1.0
        #bool_slice[x1:x2,y1:y2]=1
        print "    done slicing"

        #self.final_array += image_slice
        self.final_array[x1:x2,y1:y2] += image_cutout

        print "    done adding image to final array"

        #finite_bool = np.isfinite(self.final_array)
        #num_infinite = np.where(finite_bool==False)[0].shape[0]

        #print "    Final array INF count and max:", num_infinite, np.max(self.final_array)

        print " "
        return 1 #sunrise_image,scalesim


    def write_success_table(self,filename):
        boolthing = np.ones_like(hdudf_h.mstar_list)
        i_fail = np.where(np.asarray(hdudf_h.image_files)=='')[0]
        print boolthing.shape, i_fail.shape
        print boolthing, i_fail
        boolthing[i_fail] = 0
        data = np.asarray([boolthing,hdudf_h.x_array,hdudf_h.y_array])
        asciitable.write(data,filename)
        return

    def place_image(self,x,y,galaxy_image,galaxy_pixsize):
        new_image = self.final_array

        return new_image


if __name__=="__main__":
    mstar_list = np.asarray([8.0])
    redshift_list = np.asarray([2.0])
    sfr_list = np.asarray([1.5])
    #instead, use observed UDF catalogs

    udf_hdulist = pyfits.open('data/udf_zbest_sedfit_jen2015.fits')
    udf_table = udf_hdulist[1].data
    udf_zbest = udf_table.field('Z_BEST')
    udf_lmstar = udf_table.field('LMSTAR_BC03')
    udf_hmag = udf_table.field('MAG_F160W')

    x_list = np.asarray([2500])
    y_list = np.asarray([6000])
    #random positions?

    fi = np.where(udf_hmag > 27.0)[0]
    fake_zs = np.asarray([udf_zbest[fi],udf_zbest[fi]]).flatten()
    fake_lmasses = 7.0 - 1.6*np.random.random(fake_zs.shape[0])
    fake_hmag = 29.0 + 4.0*np.random.random(fake_zs.shape[0])

    udf_zbest = np.append(udf_zbest,fake_zs)
    udf_lmstar = np.append(udf_lmstar,fake_lmasses)
    udf_hmag = np.append(udf_hmag,fake_hmag)

    #Npix = 27800.0/2.0
    Npix = 10000.0 #16880 w/ 2.78 arcmin gives 10mas pixels
    Narcmin = 1.0
    Narcsec = Narcmin*60.0

    #Npix_hst = 27800.0/8.0 #ish
    Npix_hst = 1200.0

    Pix_arcsec = Narcsec/Npix
    Pix_arcsec_hst = Narcsec/Npix_hst

    print "Modeling image with pixel scale (arcsec): ", Pix_arcsec

    blank_array = np.float32(np.zeros(shape=(Npix,Npix)))
    print blank_array.shape
    blank_array_hst = np.float32(np.zeros(shape=(Npix_hst,Npix_hst)))
    
    sim_catalog_file = '/Users/gsnyder/Documents/Projects/HydroART_Morphology/Hyades_Data/juxtaposicion-catalog-Nov18_2013/data/sim'
    simdata = asciitable.read(sim_catalog_file,data_start=1)
    #print simdata

    rf = ['F850LP','F606W','F160W','F775W','F125W']

    hdudf_h = mock_hdudf(Npix,Pix_arcsec,blank_array,'F160W',simdata,Narcmin,1.60,28.0,0.025,req_filters=rf)
    hdudf_j = mock_hdudf(Npix,Pix_arcsec,blank_array,'F125W',simdata,Narcmin,1.25,28.0,0.022,req_filters=rf)
    hdudf_z = mock_hdudf(Npix,Pix_arcsec,blank_array,'F850LP',simdata,Narcmin,0.90,28.0,0.015,req_filters=rf)
    hdudf_i = mock_hdudf(Npix,Pix_arcsec,blank_array,'F775W',simdata,Narcmin,0.75,28.0,0.014,req_filters=rf)
    hdudf_v = mock_hdudf(Npix,Pix_arcsec,blank_array,'F606W',simdata,Narcmin,0.60,28.0,0.012,req_filters=rf)

    hudf_h = mock_hdudf(Npix_hst,Pix_arcsec_hst,blank_array_hst,'F160W',simdata,Narcmin,1.60,28.0,0.20,req_filters=rf)
    hudf_z = mock_hdudf(Npix_hst,Pix_arcsec_hst,blank_array_hst,'F850LP',simdata,Narcmin,0.90,28.0,0.15,req_filters=rf)
    hudf_v = mock_hdudf(Npix_hst,Pix_arcsec_hst,blank_array_hst,'F606W',simdata,Narcmin,0.60,28.0,0.12,req_filters=rf)

    udf_success = np.int32(np.zeros_like(udf_lmstar))

    for i,z in enumerate(udf_zbest):
        if i > 50000:
            continue
        if i % 4 != 0:
            continue

        #i = udf_zbest.shape[0] - i - 1

        result = hdudf_h.find_image(udf_lmstar[i],z,0.0,i,-1,-1,udf_hmag[i])
        udf_success[i] = result
        result = hdudf_j.find_image(udf_lmstar[i],z,0.0,i,-1,-1,udf_hmag[i])
        result = hdudf_z.find_image(udf_lmstar[i],z,0.0,i,-1,-1,udf_hmag[i])
        result = hdudf_i.find_image(udf_lmstar[i],z,0.0,i,-1,-1,udf_hmag[i])
        result = hdudf_v.find_image(udf_lmstar[i],z,0.0,i,-1,-1,udf_hmag[i])

        hudf_h.find_image(udf_lmstar[i],z,0.0,i,-1,-1,udf_hmag[i])
        hudf_z.find_image(udf_lmstar[i],z,0.0,i,-1,-1,udf_hmag[i])
        hudf_v.find_image(udf_lmstar[i],z,0.0,i,-1,-1,udf_hmag[i])

        #NOTE NOW RETURNS IN Jy/pix!!

    i_fail = np.where(np.asarray(hdudf_h.image_files)=='')[0]
    print "Numfail: ", i_fail.shape

    print udf_lmstar[0:100]
    print udf_success[0:100]

    successes = {'udf_z':udf_zbest, 'udf_lmstar':udf_lmstar, 'mockudf_success':udf_success}
    ascii.write(successes,'hdudf_success_list.txt')

    #exit()

    #im,pix_pc = hdudf_h.find_image(mstar_list[0],redshift_list[0],sfr_list[0],1,x_list[0],y_list[0])
    #im2 = hdudf_h.modify_and_place(im,x_list[0],y_list[0],redshift_list[0])

    print hdudf_h.image_files
    print hdudf_h.N_inf



    #WANT ability to know which UDF entries were successful -- save image files.  Pickle?  FITS table?

    print np.max(hdudf_h.final_array)
    new_float = np.float32(hdudf_h.final_array)
    print np.max(new_float)
    new_bool = np.isfinite(new_float)
    print np.where(new_bool==False)[0].shape[0]



    primhdu = pyfits.PrimaryHDU(new_float) ; primhdu.header['IMUNIT']=('Jy/pix') ; primhdu.header['PIXSCALE']=(Pix_arcsec, 'arcsec')
    hdulist = pyfits.HDUList([primhdu])
    hdulist.writeto('hdudf_6mas_F160W_Jy.fits',clobber=True)

    primhdu = pyfits.PrimaryHDU(np.float32(hdudf_j.final_array)) ; primhdu.header['IMUNIT']=('Jy/pix') ; primhdu.header['PIXSCALE']=(Pix_arcsec, 'arcsec')
    hdulist = pyfits.HDUList([primhdu])
    hdulist.writeto('hdudf_6mas_F125W_Jy.fits',clobber=True)

    primhdu = pyfits.PrimaryHDU(np.float32(hdudf_z.final_array)) ; primhdu.header['IMUNIT']=('Jy/pix') ; primhdu.header['PIXSCALE']=(Pix_arcsec, 'arcsec')
    hdulist = pyfits.HDUList([primhdu])
    hdulist.writeto('hdudf_6mas_F850LP_Jy.fits',clobber=True)

    primhdu = pyfits.PrimaryHDU(np.float32(hdudf_i.final_array)) ; primhdu.header['IMUNIT']=('Jy/pix') ; primhdu.header['PIXSCALE']=(Pix_arcsec, 'arcsec')
    hdulist = pyfits.HDUList([primhdu])
    hdulist.writeto('hdudf_6mas_F775W_Jy.fits',clobber=True)

    primhdu = pyfits.PrimaryHDU(np.float32(hdudf_v.final_array)) ; primhdu.header['IMUNIT']=('Jy/pix') ; primhdu.header['PIXSCALE']=(Pix_arcsec, 'arcsec')
    hdulist = pyfits.HDUList([primhdu])
    hdulist.writeto('hdudf_6mas_F606W_Jy.fits',clobber=True)



    primhdu = pyfits.PrimaryHDU(np.float32(hudf_h.final_array)) ; primhdu.header['IMUNIT']=('nJy/Sr') ; primhdu.header['PIXSCALE']=(Pix_arcsec_hst, 'arcsec')
    hdulist = pyfits.HDUList([primhdu])
    hdulist.writeto('hudf_F160W_Jy.fits',clobber=True)
    primhdu = pyfits.PrimaryHDU(np.float32(hudf_z.final_array)) ; primhdu.header['IMUNIT']=('nJy/Sr') ; primhdu.header['PIXSCALE']=(Pix_arcsec_hst, 'arcsec')
    hdulist = pyfits.HDUList([primhdu])
    hdulist.writeto('hudf_F850LP_Jy.fits',clobber=True)
    primhdu = pyfits.PrimaryHDU(np.float32(hudf_v.final_array)) ; primhdu.header['IMUNIT']=('nJy/Sr') ; primhdu.header['PIXSCALE']=(Pix_arcsec_hst, 'arcsec')
    hdulist = pyfits.HDUList([primhdu])
    hdulist.writeto('hudf_F606W_Jy.fits',clobber=True)


    #hdudf_h.write_success_table('F160W_successes.txt')



    #b=hdudf_v.final_array
    #g=hdudf_z.final_array
    #r=hdudf_h.final_array
    #res = render_hdudf(b,g,r,'HDUDF_v1.pdf',pixel_arcsec=Pix_arcsec)

