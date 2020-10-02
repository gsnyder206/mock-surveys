import numpy as np
import astropy
from astropy.io import fits
from astropy.io import ascii
from astropy import wcs
import tng_api_utils as tau
import initialize_mock_fits as imf
import os
from tqdm import tqdm
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy import units as u
import time


lc_colnames=['Snapshot number', 'Subhalo index', 'RA degree', 'DEC degree',
              'RA true z', 'DEC true z', 'RA inferred z', 'DEC inferred z',
              'True z', 'Inferred z', 'Peculiar z', 'True scale',
              'Comoving X', 'Comoving Y', 'Comoving Z',
              'True angular distance', 'Inferred angular distance',
              'Snapshot z', 'Geometric z', 'Lightcone number',
              'Stellar mass w2sr', 'Total gas mass w2sr', 'Total subhalo mass',
              'Total BH mass w2sr', 'Total baryon mass w2sr', 'SFR w2sr',
              'Total BH accretion rate', 'Camera X', 'Camera Y', 'Camera Z',
              'Intrinsic g mag', 'Intrinsic r mag', 'Intrinsic i mag',
              'Intrinsic z mag', 'Galaxy motion X', 'Galaxy motion Y',
              'Galaxy motion Z/Peculiar', 'Cosmological expansion',
              'Apparent total gmag']

pTdict={'mstar':'stars','mgas':'gas','sfr':'gas','zgas':'gas','zstar':'stars',
        'halpha':'gas','hbeta':'gas','o3':'gas','n2':'gas','age':'stars','mdm':'dm'}
pFdict={'mstar':'mass','mgas':'mass','sfr':'sfr','zgas':'metal','zstar':'metal',
        'halpha':'sb_H-alpha','hbeta':'sb_H-beta','o3':'sb_O--3-5006.84A','n2':'sb_N--2-6583.45A',
        'age':'stellar_age','mdm':'mass'}
unitdict={'mstar':'Msun','mgas':'Msun','sfr':'Msun/year','zgas':'MZ/Mtot','zstar':'MZ/Mtot',
        'halpha':'photons/s/cm^2/arcsec^2','hbeta':'photons/s/cm^2/arcsec^2','o3':'photons/s/cm^2/arcsec^2','n2':'photons/s/cm^2/arcsec^2',
        'age':'Gyr','mdm':'Msun'}


def populate_hydro_source_only(hydro_cutout,cutout_size,scale_arcsec,
                                simname,snapnum,subhalo_id,pT_use,pF_use,key_use):

    axes=['0,1','0,2','1,2','1,0','2,0','2,1']

    if cutout_size > 2000:
        npix=2000
        reshape=True
    else:
        npix=cutout_size
        reshape=False

    size_arcmin=np.float32(cutout_size)*scale_arcsec/60.0

    wcs_header=hydro_cutout.wcs.to_header()



    use_hdu = tau.get_subhalo_mockdata_as_fits(sim=simname,
                                                snap=np.int32(snapnum),
                                                sfid=np.int64(subhalo_id),
                                                partType=pT_use,
                                                partField=pF_use,
                                                size=size_arcmin,
                                                nPixels=npix,
                                                axes=axes[0],
                                                existingheader=wcs_header)
    use_hdu.header['KEY']=key_use


    #manage NANs here???

    use_hdu.data[np.isnan(use_hdu.data)]=0.0

    #manage weird-unit things here like Mz/Mtot and Age

    if reshape is True:
        total_flux_njy=np.sum(use_hdu.data)
        new_data=imresize(use_hdu.data,(cutout_size,cutout_size),interp='bicubic')
        new_data = new_data*total_flux_njy/np.sum(new_data)
        use_hdu = fits.ImageHDU(new_data,header=use_hdu.header)

    input_data=hydro_cutout.data
    input_data += use_hdu.data


    #output_data = input_data*0.0 + use_hdu.data
    #input_data += output_data

    n_arcmin=use_hdu.header['N_ARCMIN']
    total_quant = np.sum(use_hdu.data)

    return n_arcmin, total_quant


def run(lcfile,filtername='wfc3_ir_f160w',outfile='test.fits',**kwargs):

    lcfname=os.path.basename(lcfile)
    simname=lcfname.split('_')[1]

    print(lcfname)
    print(simname)


    #add quantity/units functionality
    if filtername in pTdict.keys():
        pT_use=pTdict[filtername]
        pF_use=pFdict[filtername]
        #bN_use=''
        key_use=filtername
    else:
        pT_use='stars'
        pF_use='stellarBandObsFrame-'+filtername
        bN_use=filtername
        key_use=filtername


    if os.path.lexists(outfile):
        print('output file exists, picking up mid-stream...')
        hdulist=fits.open(outfile)
        hydro_hdu=hdulist[filtername]
        orig_lctable=hdulist['InputData']
        image_catalog=hdulist['OutputData']
        start_i=hydro_hdu.header['NEXTI']
    else:
        hydro_hdu = imf.blank_image(header_only=False,**kwargs)
        hydro_hdu.header['EXTNAME']=filtername
        hydro_hdu.header['FILTER']=filtername
        hydro_hdu.header['SIMNAME']=simname
        hydro_hdu.header['LCFNAME']=(lcfname,'lightcone file path')

        if filtername in unitdict.keys():
            this_bunit = unitdict[filtername]
            hydro_hdu.header['BUNIT']=this_bunit
        else:
            hydro_hdu.header['BUNIT']='nanoJanskies'



        start_i=0
        lcdata=ascii.read(lcfile,names=lc_colnames)
        #print(lcdata.info())

        orig_lctable=fits.BinTableHDU(lcdata)
        orig_lctable.header['EXTNAME']='InputData'

        orig_cols=orig_lctable.columns
        new_cols=fits.ColDefs([
            fits.Column(name='image_success',format='K',array=np.zeros_like(orig_lctable.data['Snapshot number'])),
            fits.Column(name='photrad_kpc',format='D',array=np.zeros_like(orig_lctable.data['Apparent total gmag'])),
            fits.Column(name='cutoutfov_kpc',format='D',array=np.zeros_like(orig_lctable.data['Apparent total gmag'])),
            fits.Column(name='cutout_size',format='K',array=np.zeros_like(orig_lctable.data['Snapshot number'])),
            fits.Column(name='n_arcmin',format='D',array=np.zeros_like(orig_lctable.data['Apparent total gmag'])),
            fits.Column(name='total_quant',format='D',array=np.zeros_like(orig_lctable.data['Apparent total gmag']))])

        image_catalog=fits.BinTableHDU.from_columns(new_cols+orig_cols)
        image_catalog.header['EXTNAME']='OutputData'

    scale_arcsec=hydro_hdu.header['PIXSCALE']
    hydro_wcs=wcs.WCS(hydro_hdu.header)

    for i,row in enumerate(tqdm(orig_lctable.data[start_i:start_i+110 ],mininterval=1,miniters=10,smoothing=0.1)):
        snapnum=row['Snapshot number']
        subhalo_id=row['Subhalo index']
        #print(snapnum,subhalo_id)

        #obtain size estimate
        try:
            subhalo_url=tau.baseUrl+simname+'/snapshots/'+str(snapnum)+'/subhalos/'+str(subhalo_id)
            r=tau.get(subhalo_url)
            photrad_kpc=r['stellarphotometricsrad']/tau.tngh
            starrad_kpc=r['subhalohalfmassradtype'][4]/tau.tngh
            #check subhalo flag here?

        except:
            print('failed to get subhalo info',snapnum,subhalo_id)
            continue

        #probably need something better here.
        cutoutfov_kpc=25*starrad_kpc


        kpc_per_arcsec=tau.tngcos.kpc_proper_per_arcmin(row['True z']).value/60.0

        #cut out at least 10x10
        cutout_size=np.max([10, np.int64( (cutoutfov_kpc/kpc_per_arcsec)/scale_arcsec)])
        cutout_pos=SkyCoord(row['RA degree']*u.deg,row['DEC degree']*u.deg)


        image_catalog.data['photrad_kpc'][start_i+i]=photrad_kpc
        image_catalog.data['cutoutfov_kpc'][start_i+i]=cutoutfov_kpc
        image_catalog.data['cutout_size'][start_i+i]=cutout_size

        #create cutout object
        try:
            hydro_cutout=Cutout2D(hydro_hdu.data,cutout_pos,cutout_size,wcs=hydro_wcs,mode='strict')
        except astropy.nddata.utils.NoOverlapError as NOE:
            print('No overlap between source and field image, skipping:', subhalo_url)
            continue
        except astropy.nddata.utils.PartialOverlapError as POE:
            print('Partial overlap between source and field image, skipping:', subhalo_url)
            continue

        #populate cutout object and add to images

        #populate_hydro_source_only(hydro_cutout,cutout_size,scale_arcsec,filtername,simname,snapnum,subhalo_id)


        try:
            n_arcmin,total_quant=populate_hydro_source_only(hydro_cutout,cutout_size,scale_arcsec,
                                                                    simname,snapnum,subhalo_id,
                                                                    pT_use,pF_use,key_use) #first try
        except:
            print(sys.exc_info())
            time.sleep(300) #likely server failure/reboot -- wait 5 minutes
            try:
                n_arcmin,total_quant=populate_hydro_source_only(hydro_cutout,cutout_size,scale_arcsec,
                                                                    simname,snapnum,subhalo_id,
                                                                    pT_use,pF_use,key_use)  #2nd try
            except Exception as EXC:
                print(sys.exc_info())
                print('Two failures to populate hydro source (likely server issue?), halting...')
                break #admit defeat

        image_catalog.data['n_arcmin'][start_i+i]=n_arcmin
        image_catalog.data['total_quant'][start_i+i]=total_quant
        image_catalog.data['image_success'][start_i+i]=1

        if i % 1000 == 0:
            print('Saving intermediate image')
            hydro_hdu.header['NEXTI']=start_i+i+1
            out_hdulist=fits.HDUList([hydro_hdu,orig_lctable,image_catalog])
            out_hdulist.writeto(outfile,overwrite=True)


    print('finished populating image')
    #ensure unit correctness
    hydro_hdu.header['NEXTI']=start_i+i+1
    out_hdulist=fits.HDUList([hydro_hdu,orig_lctable,image_catalog])
    out_hdulist.writeto(outfile,overwrite=True)

    return hydro_hdu
