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
import datetime

import illustris_api_utils as iau


cat_dict={'snapshot' 	      :'Illustris snapshot number',
          'SubfindID'		      :'Illustris subhalo Subfind ID number',
          'ra_deg'		      :'arbitrary right ascension in degrees',
          'dec_deg'		      :'arbitrary declination in degrees',
          'ra_kpc'		      :'arbitrary x offset in physical kpc at true redshift',
          'dec_kpc'		      :'arbitrary y offset in physical kpc at true redshift',
          'ra_kpc_inferred'	      :'same as previous but for inferred redshift',
          'dec_kpc_inferred'	      :'--',
          'true_z'		      :'true cosmological redshift',
          'inferred_z'		      :'inferred cosmological redshift (true z + peculiar v/c)           		       ',
          'peculiar_z'		      :'peculiar v/c',
          'true_kpc_per_arcsec'      :'kpc per arcsec at true redshift',
          'X_cmpc'		      :'galaxy X position in observer frame [comoving Mpc]',
          'Y_cmpc'		      :'galaxy Y position in observer frame [comoving Mpc]',
          'Z_cmpc'		      :'galaxy Z position in observer frame [comoving Mpc]',
          'ADD_cmpc'		      :'angular diameter distance at true redshift [comoving Mpc]',
          'ADD_cmpc_inferred'	      :'angular diameter distance at inferred redshift [comoving Mpc]',
          'snapshot_z'		      :'redshift of Illustris snapshot from which this subhalo is drawn',
          'geometric_z'	      :'estimated redshift at the comoving distance in the center of the lightcone path through this snapshot',
          'cylinder_number'	      :'index of snapshot repetition for this subhalo',
          'mstar_msun_rad'	      :'stellar mass in 2X stellar half-mass radius [Msun]',
          'mgas_msun_rad'	      :'gas mass in 2X stellar half-mass radius [Msun]',
          'subhalo_mass_msun'	      :'total mass in subhalo [Msun]',
          'bhmass_msun_rad'	      :'supermassive black hole mass in 2X stellar half-mass radius [Msun]',
          'mbary_msun_rad'	      :'total baryon mass in 2X stellar half-mass radius [Msun]',
          'sfr_msunperyr_rad'	      :'SFR in 2X stellar half-mass radius [Msun]',
          'bhrate_code'	      :'supermassive black hole accretion rate in code units (see data release docs)',
          'camX_mpc'		      :'galaxy X position in camera frame [physical Mpc]',
          'camY_mpc'		      :'galaxy Y position in camera frame [physical Mpc]',
          'camZ_mpc'		      :'galaxy Z position in camera frame [physical Mpc]',
          'g_AB_absmag'	      :'rest-frame SDSS-g AB absolute magnitude',
          'r_AB_absmag'	      :'rest-frame SDSS-r AB absolute magnitude',
          'i_AB_absmag'	      :'rest-frame SDSS-i AB absolute magnitude',
          'z_AB_absmag'	      :'rest-frame SDSS-z AB absolute magnitude',
          'v_kms_camX'		      :'transverse x velocity [km/s]',
          'v_kms_camY'		      :'transverse y velocity [km/s]',
          'v_kms_camZ'		      :'line-of-sight (peculiar) velocity [km/s]',
          'v_kms_hubble'	      :'true cosmological expansion velocity [km/s]',
          'g_AB_appmag'	      :'rest-frame SDSS-g AB apparent magnitude',
          'sim'		      :'simulation name (e.g., Illustris-1)',
          'snap'		      :'Illustris snapshot number (same as snapshot)',
          'sfid'		      :'Illustris subhalo Subfind ID number (same as SubfindID)',
          'z'			      :'true cosmological redshift (same as true_z)',
          'RA'			      :'same as ra_deg',
          'DEC'		      :'same as dec_deg',
          'origin_i'		      :'Internal variable describing the image array location to place the origin of this galaxy image [Internal-- does not describe final image]',
          'origin_j'		      :'Internal variable describing the image array location to place the origin of this galaxy image [Internal-- does not describe final image]',
          'pos_i'		      :'position of this subhalo in array rows  [it appears that (pos_j,pos_i) is the correct placement in Python terms] [Internal-- does not describe final image]',
          'pos_j'		      :'position of this subhalo in array columns  [it appears that (pos_j,pos_i) is the correct placement in Python terms] [Internal-- does not describe final image]',
          'pixsize_arcsec'	      :'pixel size of image [Arcsec-- Internal--does not describe final image]',
          'final_fov_arcsec'	      :'final field-of-view of image [Arcsec]',
          'full_npix'		      :'total number of pixels across image [Internal-- does not describe final image]',
          'this_npix'		      :'number of pixels in this galaxy sub-image [Internal-- does not describe final image]',
          'this_fov_kpc'	      :'field-of-view of this galaxy image in physical kpc',
          'halfmassrad_factor'	      :'factor describing this galaxy field of view in units of its stellar half-mass radius',
          'nrays'		      :'number of rays used to render this galaxy image',
          'run_dir'		      :'internal run directory for locating this galaxy sub-image',
          'success'		      :'whether or not this galaxy is rendered in the final image (should all be True; Falses were dropped from this table before saving)',
          'new_i'		      :'position of galaxy in image pixel units [final version-- use this one]  [it appears that (pos_j,pos_i) is the correct placement in Python terms]',
          'new_j'		      :'position of galaxy in image pixel units [final version-- use this one]  [it appears that (pos_j,pos_i) is the correct placement in Python terms]',
          'AB_absmag'                 :'total absolute magnitude (AB system) of galaxy in this filter'}


if __name__=="__main__":
    files = np.asarray(glob.glob('*.fits'))

    #simdata=iau.get('http://www.illustris-project.org/api/Illustris-1')

    api_web='http://www.illustris-project.org/data/docs/api/'
    
    for fn in files:
        hdulist=pyfits.open(fn)

        #update main header
        hdulist['IMAGE_NOPSF'].header['DOI']='https://doi.org/10.1093/mnras/stx487'
        hdulist['IMAGE_NOPSF'].header['AUTHOR']='Gregory F. Snyder'
        hdulist['IMAGE_NOPSF'].header['PAPER']='Snyder et al. 2017, MNRAS, 468, 207'
        hdulist['IMAGE_NOPSF'].header['DATE']=datetime.datetime.now().isoformat()

        try:
            oldfil=hdulist['IMAGE_NOPSF'].header['ALLFIL']
        except:
            oldfil=hdulist['IMAGE_NOPSF'].header['FILTER']

        split1=oldfil.split('-')
        telescope_sub=split1[0]
        instrument_sub=split1[1].split('_')[0]
        filter_sub=split1[1].split('_')[1]
        print(telescope_sub, instrument_sub, filter_sub)

        hdulist['IMAGE_NOPSF'].header['MISSION']=(telescope_sub.upper(),'Mission/telescope')
        hdulist['IMAGE_NOPSF'].header['INSTR']=(instrument_sub.upper(),'Instrument')
        hdulist['IMAGE_NOPSF'].header['TELESCOPE']=(telescope_sub.upper(),'Mission/telescope')
        hdulist['IMAGE_NOPSF'].header['INSTRUMENT']=(instrument_sub.upper(),'Instrument')
        hdulist['IMAGE_NOPSF'].header['FILTER']=(filter_sub.upper(),'filter')
        hdulist['IMAGE_NOPSF'].header['ALLFIL']=oldfil
        hdulist['IMAGE_NOPSF'].header['SIM_NAME']='Illustris-1'
        hdulist['IMAGE_NOPSF'].header['SIM_DATA']='http://www.illustris-project.org'
        hdulist['IMAGE_NOPSF'].header['IMTYPE']=('Survey','type of image source')
        
        hdulist['IMAGE_NOPSF'].header['REDSHIFT']=('0.5-20','Redshift of object or survey')
        pix_nopsf = hdulist['IMAGE_NOPSF'].header['PIXSIZE']
        units_nopsf=hdulist['IMAGE_NOPSF'].header['UNIT']

        try:
            #add info to all image headers
            hdulist['IMAGE_PSF'].header['REDSHIFT']=('0.5-20','Redshift of object or survey')
            hdulist['IMAGE_PSF'].header['PIXSIZE']=(pix_nopsf,'arcsec')
            hdulist['IMAGE_PSF'].header['UNIT']=(units_nopsf,'per pixel')
            hdulist['MODELPSF'].header['PIXSIZE']=(pix_nopsf,'arcsec')
        except:
            pass
        
            
        
        #add sim parameter header/table
        sim_dict={'url':'http://www.illustris-project.org/api/Illustris-1',
                  'apidoc':'http://www.illustris-project.org/data/docs/api/'}
        sim_table=astropy.table.Table([sim_dict])
        sim_hdu=pyfits.table_to_hdu(sim_table)

        sim_hdu.header['EXTNAME']='SimulationAssumptions'
        sim_hdu.header['APIDOC']=('illustris-project.org/data/docs/api/','Illustris Data API Docs')
        sim_hdu.header['URL']=('illustris-project.org/api/Illustris-1','Simulation Parameters')
        


        #add mock data header/table
        mock_hdu=pyfits.ImageHDU()
        mock_hdu.header['EXTNAME']='MockDataAssumptions'
        mock_hdu.header['CODE']='Sunrise'
        mock_hdu.header['SMODEL']='Starburst99'
        mock_hdu.header['IMF']='Kroupa'
        mock_hdu.header['Zs']=('Multiple','stellar metallicities')
        mock_hdu.header['DUST']='None'
        mock_hdu.header['SMOOTH']=('NGB64','see Torrey et al. 2015')
        
        

        #add info to catalog table header ?
        catdoc_table=astropy.table.Table([cat_dict])
        catdoc_hdu=pyfits.table_to_hdu(catdoc_table)
        catdoc_hdu.header['EXTNAME']='CatalogDocumentation'
        hdulist.append(catdoc_hdu)

        try:
            new_hdulist=pyfits.HDUList([hdulist['IMAGE_NOPSF'],sim_hdu,mock_hdu,hdulist['IMAGE_PSF'],hdulist['MODELPSF'],hdulist['Catalog'],catdoc_hdu])
        except:
            new_hdulist=pyfits.HDUList([hdulist['IMAGE_NOPSF'],sim_hdu,mock_hdu,hdulist['Catalog'],catdoc_hdu])

        
        new_hdulist.writeto(hdulist.filename(),overwrite=True)

         
