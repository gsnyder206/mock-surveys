import create_mock_hydro_image as cmhi
import sys

#'wfc3_ir_f160w'
#jwst_f200w
#wfc_acs_f814w

fdict={
'F435W':'wfc_acs_f435w',
'F606W':'wfc_acs_f606w',
'F775W':'wfc_acs_f775w',
'F814W':'wfc_acs_f814w',
'F850LP':'wfc_acs_f850lp',

'F105W':'wfc3_ir_f105w',
'F125W':'wfc3_ir_f125w',
'F140W':'wfc3_ir_f140w',
'F160W':'wfc3_ir_f160w',

'F070W':'jwst_f070w',  #non-nyquist
'F090W':'jwst_f090w',  #non-nyquist
'F115W':'jwst_f115w',  #non-nyquist
'F150W':'jwst_f150w',  #non-nyquist
'F200W':'jwst_f200w',
'F277W':'jwst_f277w',
'F356W':'jwst_f356w',
'F444W':'jwst_f444w'
}

pdict={
'F435W':0.03,
'F606W':0.03,
'F775W':0.03,
'F814W':0.03,
'F850LP':0.03,

'F105W':0.06,
'F125W':0.06,
'F140W':0.06,
'F160W':0.06,

'F070W':0.03,  #non-nyquist
'F090W':0.03,  #non-nyquist
'F115W':0.03,  #non-nyquist
'F150W':0.03,  #non-nyquist
'F200W':0.03,
'F277W':0.03,
'F356W':0.03,
'F444W':0.03
}

# Lightcone_TNG100-1_mag30_7_6_xyz.txt

if len(sys.argv) != 2:
    print('usage: python runscript FKEY')
    print('Available FKEY values:  ')
    print(fdict.keys())
    exit()

fkey=sys.argv[1]
assert fkey in fdict.keys()
assert fkey in pdict.keys()


hh=cmhi.run('priority_lightcones/Lightcone_TNG100-1_mag30_7_6_xyz.txt',
            scale_arcsec=pdict[fkey],w_deg=0.2,h_deg=0.2,
            filtername=fdict[fkey],
            outfile='outputs/TNG100-1_xyz_'+fkey+'.fits')


#hh=cmhi.run('priority_lightcones/testlc_TNG100-1_small.txt',
#            scale_arcsec=pdict[fkey],w_deg=0.2,h_deg=0.2,
#            filtername=fdict[fkey],
#            outfile='testlc_'+fkey+'.fits')
