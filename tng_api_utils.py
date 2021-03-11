import os
import sys
import numpy as np
import glob
import shutil
import tarfile
import string
import astropy
import astropy.cosmology
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import astropy.table
import re
import requests
import io
import h5py
import load_api_key as lak
from astropy.stats import gaussian_fwhm_to_sigma
import scipy.ndimage
import scipy as sp
import photutils
from photutils import aperture_photometry
from photutils import CircularAperture


ilh = 0.704
illcos = astropy.cosmology.FlatLambdaCDM(H0=70.4,Om0=0.2726,Ob0=0.0456)

tngh = 0.6774
tngcos = astropy.cosmology.FlatLambdaCDM(H0=100.0*tngh,Om0=0.3089,Ob0=0.0486)

illustris_defaultparams={'stars':'ParticleIDs,Coordinates,Velocities,GFM_StellarFormationTime,GFM_Metallicity,GFM_InitialMass,Masses','gas':'ParticleIDs,Coordinates,Density,ElectronAbundance,Masses,Velocities,Volume,SubfindDensity,Potential,InternalEnergy,StarFormationRate,GFM_Metallicity,GFM_AGNRadiation,GFM_WindDMVelDisp,GFM_CoolingRate,NeutralHydrogenAbundance,SmoothingLength,SubfindHsml,SubfindVelDisp,NumTracers','dm':'Coordinates,Velocities,Potential,ParticleIDs','bhs':'all'}


baseUrl = 'http://www.tng-project.org/api/'

headers = {"api-key":lak.aks}


def get(path, params=None,savepath=None, stream=False):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)
    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers and stream is False:
        file_basename = r.headers['content-disposition'].split("filename=")[1]

        if savepath is not None:
            filename = os.path.join(savepath,file_basename)
        else:
            filename = file_basename

        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string
    elif 'content-disposition' in r.headers and stream is True:
        file_basename = r.headers['content-disposition'].split("filename=")[1]

        return r # return the object

    return r


#takes subhalo and partField as inputs
#for ObsFrame, returns apparent AB magnitudes
#for rest frame, returns absolute AB magnitudes

def get_subhalo_magnitude(sim='TNG100-1',
                            snap=50,
                            sfid=0,
                            partField='stellarBandObsFrame-wfc3_ir_f160w',
                            axes='0,1'):


    partType='stars'
    size=0.30
    nPixels=1

    url=baseUrl+'/'+sim+'/snapshots/'+str(snap)+'/subhalos/'+str(sfid)+\
             '/vis.hdf5?partType='+partType+'&partField='+partField+'&size='+str(size)+\
             '&sizeType=arcmin&nPixels='+str(nPixels)+'&method=sphMap_subhalo'+'&axes='+axes



    r=get(url,stream=True)

    output=io.BytesIO()
    try:
        output.write(r.content)
        ho=h5py.File(output,mode='r')
    except OSError as OE:
        print(r)
        print(url)
        raise

    #print(ho.keys())
    #print(ho['Header'])
    data=ho['grid'][()] #dataset should be nPixels x nPixels

    pixsize_arcsec = 60.0*(size/nPixels)

    if partField.find('stellarBandObsFrame-')==0:
        in_units='mag/arcsec^2' #(confirmed)
        #convert to nanoJanskies
        #size is fixed to be in arcmin
        flux_njy=(pixsize_arcsec**2)*1.0e9*3631.0*np.sum(10.0**(-0.4*data))  #nJy
        magnitude= -2.5*np.log10(flux_njy*(1.0e-9)/3631.0)

    elif partField.find('stellarBand-')==0:
        in_units='ABS AB mag'  #(true?)
        flux= np.sum(10.0**(-0.4*data))
        magnitude= -2.5*np.log10(flux)
    else:
        assert(False)

    return magnitude



def get_subhalo_mockdata_as_fits(sim='TNG100-1',
                                 snap=50,
                                 sfid=0,
                                 partType='stars',
                                 partField='stellarBandObsFrame-',
                                 size=0.1,
                                 nPixels=150,
                                 axes='0,1',
                                 existingheader=None):


    if nPixels > 2000:
        return None
    #+'&axes='+axes
    url=baseUrl+'/'+sim+'/snapshots/'+str(snap)+'/subhalos/'+str(sfid)+\
             '/vis.hdf5?partType='+partType+'&partField='+partField+'&size='+str(size)+'&sizeType=arcmin&nPixels='+str(nPixels)+'&method=sphMap_subhalo'+'&axes='+axes



    r=get(url,stream=True)

    output=io.BytesIO()
    try:
        output.write(r.content)
        ho=h5py.File(output,mode='r')
    except OSError as OE:
        print(r)
        print(url)
        raise

    #print(ho.keys())
    #print(ho['Header'])
    data=ho['grid'][()] #dataset should be nPixels x nPixels

    pixsize_arcsec = 60.0*(size/nPixels)

    if partField.find('sb_')==0:
        in_units='photons/s/cm^2/arcsec^2'
        out_data = data*1.0
        out_units=in_units
    elif partField.find('stellarBandObsFrame-')==0:
        in_units='mag/arcsec^2' #(confirmed)
        #convert to nanoJanskies
        #size is fixed to be in arcmin
        flux_njy=(pixsize_arcsec**2)*1.0e9*3631.0*10.0**(-0.4*data)  #nJy
        out_data = flux_njy*1.0
        out_units='nanoJanskies'
    elif partField.find('stellarBand-')==0:
        in_units='ABS AB mag'  #(true?)
        out_data = data*1.0
        out_units=in_units
    else:
        #need to confirm other field's units before allowing them here
        #MASS, SFR, METALLICITY MAPS OMG
        #need to figure out what to do about NANs in surface brightness maps?  also metal mass densities?


        inunitdict={'mass':'log Msun','sfr':'log Msun/year','metal':'log MZ/Mtot',
                    'stellar_age':'Gyr'}
        outunitdict={'mass':'Msun','sfr':'Msun/year','metal':'MZ/Mtot',
                    'stellar_age':'Gyr'}

        if partField in inunitdict.keys():
            in_units=inunitdict[partField]
            out_units=outunitdict[partField]
            if partField=='stellar_age':
                out_data=1.0*data #keep in Gyr
            else:
                out_data=10.0**data #de-log
        else:
            return None


    if existingheader is not None:
        fits_hdu = fits.ImageHDU(out_data,header=existingheader)
    else:
        fits_hdu = fits.ImageHDU(out_data)

    fits_hdu.header['EXTNAME']='MockData'
    fits_hdu.header['SIM']=sim
    fits_hdu.header['SNAP']=snap
    fits_hdu.header['SFID']=sfid
    fits_hdu.header['pType']=partType
    fits_hdu.header['pField']=partField
    #fits_hdu.header['bandName']=bandName
    fits_hdu.header['N_arcmin']=size
    fits_hdu.header['nPixels']=nPixels
    fits_hdu.header['axes']=axes
    fits_hdu.header['pixscale']=(pixsize_arcsec,'arcsec')
    fits_hdu.header['origunit']=(in_units,'downloaded units')
    fits_hdu.header['BUNIT']=(out_units,'final units')


    return fits_hdu


#takes as input the results of get_subhalo_mockdata_as_fits
#returns the new HDU with PSF applied

def convolve_with_fwhm(in_hdu, fwhm_arcsec=0.10):

    #load image data and metadata
    image_in=in_hdu.data
    header_in=in_hdu.header

    #Choi et al. 2018
    #this_cosmology=astropy.cosmology.FlatLambdaCDM(0.72*100.0,0.26,Ob0=0.044)
    #kpc_per_arcsec = this_cosmology.kpc_proper_per_arcmin(header_in['REDSHIFT']).value/60.0

    #calculate PSF width in pixel units
    #pixel_size_arcsec=header_in['PIX_KPC']/kpc_per_arcsec
    pixel_size_arcsec=header_in['pixscale']

    sigma_arcsec=fwhm_arcsec*gaussian_fwhm_to_sigma
    sigma_pixels=sigma_arcsec/pixel_size_arcsec

    image_out=sp.ndimage.filters.gaussian_filter(image_in,sigma_pixels,mode='nearest')
    hdu_out = fits.ImageHDU(image_out,header=header_in)
    hdu_out.header['FWHMPIX']=(sigma_pixels/gaussian_fwhm_to_sigma,'pixels')
    hdu_out.header['FWHM']=(fwhm_arcsec,'arcsec')
    hdu_out.header['SIGMAPIX']=(sigma_pixels,'pixels')
    hdu_out.header['SIGMA']=(sigma_arcsec,'arcsec')
    hdu_out.header['EXTNAME']='MockData_PSF'
    hdu_out.header['PIXSIZE']=(pixel_size_arcsec,'arcsec')

    return hdu_out


#takes as input the results of convolve_with_fwhm
#returns the new HDU with basic noise model applied
def add_simple_noise_extractedsn(in_hdu,radius_arcsec=0.5,extractedsn=300):


    image_in=in_hdu.data
    header_in=in_hdu.header

    #this is the approximate equation for sb magnitude input limit
    #sigma_njy=(2.0**(-0.5))*((1.0e9)*(3631.0/5.0)*10.0**(-0.4*sb_maglim))*header_in['PIXSIZE']*(3.0*header_in['FWHM'])

    #get flux in aperture
    #can we guarantee that the galaxy is in the middle?
    npix=image_in.shape[0]
    ci=np.float32(npix)/2

    radius_pixels = radius_arcsec/in_hdu.header['PIXSIZE']

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
