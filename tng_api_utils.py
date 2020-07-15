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
        ho=h5py.File(output)
    except OSError as OE:
        print(r)
        print(url)
        raise

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
