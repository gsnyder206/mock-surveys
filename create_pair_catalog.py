
import pandas as pd
import astropy
import astropy.io.fits as fits
import astropy.cosmology

from astropy.coordinates import SkyCoord

from astropy.table import Table
import astropy.io.ascii as ascii

import os
import numpy as np

ilh = 0.704
illcos = astropy.cosmology.FlatLambdaCDM(H0=70.4,Om0=0.2726,Ob0=0.0456)

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


if __name__=="__main__":

    #these are created with gfs_sublink_utils.py:
    
    #primaries files:
    ##iz snap sfid mstar mhalo mbary delta hsml sfr
    #pairs files:
    ##iz snap sfid mstar mhalo mbary Nc sec_iz sec_snap sec_sfid sec_mstar sec_mhalo sec_mbary sec_drad_arcsec drad_kpc pairMerges mergersnap mergerage time_until_merger mstar_ratio_now mstar_ratio_tmax mhalo_ratio_tmax mhalo_ratio_now snap_tmax merger_mass merger_mstar ecc rperi delta hsml rnow vnow rapo b pri_sfr sec_sfr delta_z delta_v gratio bratio


    #these were provided by Jen:
    
    # new/Illustris-1_RADEC_hudfwide_75Mpc_7_6_xyz_corners.txt
    # Mstar_primary=     3.1623E+10 dr_max= 142.857 kpc Mstar_c/Mstar_p = ******* delta iz/1+iz =    0.02
    # primary snapid subhaloid iz tz mstar mhalo mbary ncompanions
    # secondary  snapid subhaloid dradius(ARCSEC-CONFIRMED)  iz tz mstar  mhalo mbary 
    

    lcf=os.path.expandvars('$HOME/Dropbox/JWST_Cycle1/JWST1_Assembly_Proposal/pairs/Illustris-1_RADEC_JWST1_75Mpc_7_6_zyx.txt')
    pairf=os.path.expandvars('$HOME/Dropbox/JWST_Cycle1/JWST1_Assembly_Proposal/pairs/Pairs_7_6_zyx.txt')

    
    lcdata=ascii.read(lcf)

    for colkey in lcfile_cols:
        col_obj=lcdata[colkey]
        col_obj.name=lcfile_cols[colkey]
    
    print(lcdata)

    masslim=10.0**(9.5)
    distkpch=100.0
    distkpc=distkpch/ilh
    deltalim=0.02

    
    pri_index=np.where(lcdata['mstar_msun_rad'] >= masslim )[0]

    coords=SkyCoord(lcdata['ra_deg'],lcdata['dec_deg'],unit='deg')


    with open(pairf,'w') as pairfo:

        pairfo.write('#\n#\n#\n#\n')
        
        for pi in pri_index:
            pri_coord=SkyCoord(lcdata['ra_deg'][pi],lcdata['dec_deg'][pi],unit='deg')
            pri_iz=lcdata['inferred_z'][pi]
            
            #print(pi,pri_iz)
            
            pri_kpc_per_deg=illcos.kpc_proper_per_arcmin(pri_iz).value*60.0
            
            distances=(coords.separation(pri_coord)).value #in deg?
            distances_kpc = pri_kpc_per_deg*distances

            distances_arcsec = distances*3600.0  #deg-> arcsec
            
            delta =np.abs(lcdata['inferred_z']-pri_iz)/(1.0 + pri_iz)
            
            allsec=np.where(np.logical_and(distances_kpc <= distkpc, delta <= deltalim ))[0]
            
            #print(pri_iz, allsec.shape)

            if allsec.shape[0] > 1:
                #one primary plus all its companions
                pairfo.write('{:10d} {:10d} {:12.6f} {:12.6f} {:12.6e} {:12.6e} {:12.6e} {:10d}\n'.format(lcdata['snapshot'][pi],lcdata['SubfindID'][pi],pri_iz,lcdata['true_z'][pi],
                                                                                                        lcdata['mstar_msun_rad'][pi],lcdata['subhalo_mass_msun'][pi],lcdata['mbary_msun_rad'][pi],allsec.shape[0]-1    ))                

                for si in allsec:
                    #don't duplicate primary
                    if si==pi:
                        continue
                    
                    pairfo.write('{:10d} {:10d} {:12.6f} {:12.6f} {:12.6f} {:12.6e} {:12.6e} {:12.6e}\n'.format(lcdata['snapshot'][si],lcdata['SubfindID'][si], distances_arcsec[si], lcdata['inferred_z'][si], lcdata['true_z'][si],
                                                           lcdata['mstar_msun_rad'][si],lcdata['subhalo_mass_msun'][si],lcdata['mbary_msun_rad'][si] ))   
