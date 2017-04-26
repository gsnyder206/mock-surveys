## Illustris Mock Deep Fields
### Created by Gregory Snyder, Space Telescope Science Institute

Files with names:
* hlsp_misty_illustris_jwst-nircam_f150w_FIELDA_11_10_v1_lightcone.fits
* hlsp_misty_sim_telescope-instrument_filter_(geometry/selection)_version_type.fits

Acknowledgements.  These are based on the Illustris Project hydrodynamical simulations:
* Lightcone image presented by: Vogelsberger et al. 2014
* Data release paper and documentation:  illustris-project.org ; Nelson et al. 2015
* Lightcone application for hydro sims described by: Snyder et al. 2017
* Mock HST and JWST machinery supported by HST#13887

FITS files contain the following HDUs:
0-"IMAGE_NOPSF":  Raw image without PSF or noise.  Header documented below.
1-"IMAGE_PSF"  :  [Optional] Same image, now convolved with model PSF from TinyTim or WebbPSF.
2-"MODELPSF"   :  [Optional] The PSF model at the same pixel scale as HDUs 0 and 1
3-"Catalog"    :  Lightcone Catalog, containing intrinsic simulation info, including galaxy ID numbers and image positions

In HDU 0-"IMAGE_NOPSF", the header contains the following useful cards, for example:
FILTER  = 'jwst-nircam_f200w'                                                   
PIXSIZE =               0.0317 / arcsec                                         
UNIT    = 'nanoJanskies'       / per pixel                                      
ABZP    =    31.40006562228223 / AB mag zeropoint                               
PHOTFNU =             2.64E-08 / Jy; approx flux[Jy] at 1 count/sec             
EXTNAME = 'IMAGE_NOPSF'  

These values, especially PHOTFNU, are highly approximate and should be confirmed by the user using an ETC before applying to any analyses.

For HDUs 1-2, I produced these only for HST and JWST filters, as I have not yet fully propagated the WFIRST filters and PSF characteristics. These might be a good test bed for "in-situ" data simulation methods.  Also it might make sense to wait for the hardware to stabilize a bit more before fully installing these assumptions in a pipeline.

For HDU 3, here is the column documentation:
'snapshot' 	       	      Illustris snapshot number
'SubfindID'		      Illustris subhalo Subfind ID number
'ra_deg'		      arbitrary right ascension in degrees
'dec_deg'		      arbitrary declination in degrees
'ra_kpc'		      arbitrary x offset in physical kpc at true redshift
'dec_kpc'		      arbitrary y offset in physical kpc at true redshift
'ra_kpc_inferred'	      same as previous but for inferred redshift
'dec_kpc_inferred'	      --
'true_z'		      true cosmological redshift
'inferred_z'		      inferred cosmological redshift (true z + peculiar v/c)           		       
'peculiar_z'		      peculiar v/c
'true_kpc_per_arcsec'	      kpc per arcsec at true redshift
'X_cmpc'		      galaxy X position in observer frame [comoving Mpc]
'Y_cmpc'		      galaxy Y position in observer frame [comoving Mpc]
'Z_cmpc'		      galaxy Z position in observer frame [comoving Mpc]
'ADD_cmpc'		      angular diameter distance at true redshift [comoving Mpc]
'ADD_cmpc_inferred'	      angular diameter distance at inferred redshift [comoving Mpc]
'snapshot_z'		      redshift of Illustris snapshot from which this subhalo is drawn
'geometric_z'		      estimated redshift at the comoving distance in the center of the lightcone path through this snapshot
'cylinder_number'	      index of snapshot repetition for this subhalo
'mstar_msun_rad'	      stellar mass in 2X stellar half-mass radius [Msun]
'mgas_msun_rad'		      gas mass in 2X stellar half-mass radius [Msun]
'subhalo_mass_msun'	      total mass in subhalo [Msun]
'bhmass_msun_rad'	      supermassive black hole mass in 2X stellar half-mass radius [Msun]
'mbary_msun_rad'	      total baryon mass in 2X stellar half-mass radius [Msun]
'sfr_msunperyr_rad'	      SFR in 2X stellar half-mass radius [Msun]
'bhrate_code'		      supermassive black hole accretion rate in code units (see data release docs)
'camX_mpc'		      galaxy X position in camera frame [physical Mpc]
'camY_mpc'		      galaxy Y position in camera frame [physical Mpc]
'camZ_mpc'		      galaxy Z position in camera frame [physical Mpc]
'g_AB_absmag'		      rest-frame SDSS-g AB absolute magnitude
'r_AB_absmag'		      rest-frame SDSS-r AB absolute magnitude
'i_AB_absmag'		      rest-frame SDSS-i AB absolute magnitude
'z_AB_absmag'		      rest-frame SDSS-z AB absolute magnitude
'v_kms_camX'		      transverse x velocity [km/s]
'v_kms_camY'		      transverse y velocity [km/s]
'v_kms_camZ'		      line-of-sight (peculiar) velocity [km/s]
'v_kms_hubble'		      true cosmological expansion velocity [km/s]
'g_AB_appmag'		      rest-frame SDSS-g AB apparent magnitude

'sim'			      simulation name (e.g., 'Illustris-1')
'snap'			      Illustris snapshot number (same as 'snapshot')
'sfid'			      Illustris subhalo Subfind ID number (same as 'SubfindID')
'z'			      true cosmological redshift (same as 'true_z')
'RA'			      same as 'ra_deg'
'DEC'			      same as 'dec_deg'
'origin_i'		      Internal variable describing the image array location to place the origin of this galaxy's image [Internal-- does not describe final image]
'origin_j'		      Internal variable describing the image array location to place the origin of this galaxy's image [Internal-- does not describe final image]
'pos_i'			      position of this subhalo in array rows  [it appears that (pos_j,pos_i) is the correct placement in Python terms] [Internal-- does not describe final image]
'pos_j'			      position of this subhalo in array columns  [it appears that (pos_j,pos_i) is the correct placement in Python terms] [Internal-- does not describe final image]
'pixsize_arcsec'	      pixel size of image [Arcsec-- Internal--does not describe final image]
'final_fov_arcsec'	      final field-of-view of image [Arcsec]
'full_npix'		      total number of pixels across image [Internal-- does not describe final image]
'this_npix'		      number of pixels in this galaxy's sub-image [Internal-- does not describe final image]
'this_fov_kpc'		      field-of-view of this galaxy's image in physical kpc
'halfmassrad_factor'	      factor describing this galaxy's field of view in units of its stellar half-mass radius
'nrays'			      number of rays used to render this galaxy's image
'run_dir'		      internal run directory for locating this galaxy's sub-image
'success'		      whether or not this galaxy is rendered in the final image (should all be True; Falses were dropped from this table before saving)
'new_i'			      position of galaxy in image pixel units [final version-- use this one]  [it appears that (pos_j,pos_i) is the correct placement in Python terms]
'new_j'			      position of galaxy in image pixel units [final version-- use this one]  [it appears that (pos_j,pos_i) is the correct placement in Python terms]
