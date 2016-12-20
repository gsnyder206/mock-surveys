import math
import string
import sys
import os
import struct
import numpy as np
import astropy.io.fits as pyfits
import glob
import subprocess
import asciitable
import scipy
import scipy.signal
import scipy.ndimage
import scipy.ndimage.filters

def add_all_broadband_flux(hdu_index=12,filter_index=11,params_index=3,bb_index=1,iq_index=5,aux_index=2,dirname='SCI_IMAGES',mod=10, startz=0.0, endz=12.0, has_aux=True, default=True):

    if not os.path.exists(dirname):
        os.makedirs(dirname)


    if has_aux==True:
        print "Note: Requires AUX HDU!!"

    if (has_aux==False) and (default==True):
        iq_index = iq_index-1 ; hdu_index=hdu_index-1 ; filter_index=filter_index-1;params_index=params_index-1


    print iq_index
    
    nums = np.asarray(['0055','0056','0057','0059','0060','0061','0062','0063','0064','0065'])

    #nums=np.asarray(['0016','0017','0018','0019','0020'])
    bbfiles = np.sort(np.array(glob.glob('broadband*.fits')))#[0:10]
    
    imind = hdu_index

    hdus = pyfits.open(bbfiles[0])

    cambase = ((hdus)[imind]).data
    print cambase.shape
    singleslice = np.zeros_like(cambase[0,:,:])
    #print singleslice.shape
    totalflux = np.zeros_like(cambase) ; fluxunit = hdus[imind].header.get('IMUNIT') ; print fluxunit ; pixscale = hdus[imind].header.get('CD1_1')
    partialflux = np.zeros_like(cambase)

    totalflux_dandanxu = np.zeros_like(cambase)
    block_dandanxu = np.zeros_like(cambase)
    totalflux_gfs = np.zeros_like(cambase)
    block_gfs = np.zeros_like(cambase)

    zweightslice = np.zeros_like(cambase)
    partialzweight = np.zeros_like(cambase)
    cambasezeros = np.zeros_like(cambase)

    totalsm = np.zeros_like(singleslice) ; partialsm = np.zeros_like(singleslice)
    totalsfr = np.zeros_like(singleslice) ; partialsfr = np.zeros_like(singleslice)
    totalstellarmetals = np.zeros_like(singleslice) ; partialstellarmetals = np.zeros_like(singleslice)
    smdata = np.zeros_like(singleslice)

    
    filter_hdu = hdus[filter_index]
    params_hdu = hdus[params_index]
    iq_hdu = hdus[iq_index]

    stellarmassunit = 'Msun'
    metalsunit = 'Msun'
    if has_aux==True:
        aux_hdu = hdus[aux_index]
        stellarmassunit = aux_hdu.header.get('MS_UNIT')[0:4]#+'(kpc^2)', we're gonna multiply by the pixel area
        metalsunit = aux_hdu.header.get('MMS_UNIT')[0:4]#+'(kpc^2)'

    Nfiles = (bbfiles.shape)[0]
    Nlambda = (iq_hdu.data.field('lambda')).shape[0]
    Nfilters = (filter_hdu.data.field('filter')).shape[0]
    
    z_array = np.ndarray(shape=(Nfiles))
    pix_array = np.ndarray(shape=(Nfiles))
    camD_array = np.ndarray(shape=(Nfiles))
    fov_array = np.ndarray(shape=(Nfiles))
    sm_array = np.ndarray(shape=(Nfiles))
    mets_array = np.ndarray(shape=(Nfiles))


    total_fov = params_hdu.header.get('FOV') # in radians
    pixel_fov_arcsec = (total_fov*3600.0*180.0/math.pi ) / (1.0*singleslice.shape[0])


    L_lambda = np.ndarray(shape=(Nfiles,Nlambda)) ; L_lambda_unit = iq_hdu.header.get('TUNIT3')
    L_lambda_ns = np.ndarray(shape=(Nfiles,Nlambda)) ; L_lambda_ns_unit = iq_hdu.header.get('TUNIT5')
    L_lambda_out = np.ndarray(shape=(Nfiles,Nlambda)) ; L_lambda_out_unit = iq_hdu.header.get('TUNIT6')
    L_lambda_abs = np.ndarray(shape=(Nfiles,Nlambda)) ; L_lambda_abs_unit = iq_hdu.header.get('TUNIT4')

    lambda_eff = np.ndarray(shape=(Nfiles,Nfilters)) ; lambda_eff_unit = filter_hdu.header.get('TUNIT2')
    ewidth_lambda = np.ndarray(shape=(Nfiles,Nfilters)) ; ewidth_lambda_unit = filter_hdu.header.get('TUNIT3')
    ewidth_nu = np.ndarray(shape=(Nfiles,Nfilters)) ; ewidth_nu_unit = filter_hdu.header.get('TUNIT4')
    L_lambda_eff = np.ndarray(shape=(Nfiles,Nfilters)) ; L_lambda_eff_unit = filter_hdu.header.get('TUNIT5')
    AB_mag = np.ndarray(shape=(Nfiles,Nfilters)) 
    L_lambda_eff_ns = np.ndarray(shape=(Nfiles,Nfilters)) ; L_lambda_eff_ns_unit = filter_hdu.header.get('TUNIT7')
    AB_mag_ns = np.ndarray(shape=(Nfiles,Nfilters)) 
    intSB = np.ndarray(shape=(Nfiles,Nfilters)) ; intSB_unit = fluxunit
    cumSB = np.ndarray(shape=(Nfiles,Nfilters)) ; cumSB_unit = fluxunit


    filter_array = filter_hdu.data.field('filter')
    col = pyfits.Column(name='filter',format='A30',array=filter_array)
    filcols = pyfits.ColDefs([col])
    filname_hdu = pyfits.new_table(filcols) ; filname_hdu.update_ext_name('FilterNames')
    
    lam_array = iq_hdu.data.field('lambda') ; lamunit=iq_hdu.header.get('TUNIT1')
    col = pyfits.Column(name='lambda',format='D',unit=lamunit,array=lam_array)
    lamcols = pyfits.ColDefs([col])
    lamhdu = pyfits.new_table(lamcols) ; lamhdu.update_ext_name('Lambda')
    #lamhdu = pyfits.BinTableHDU(lamtab)



    extdata = asciitable.read('/Users/gsnyder/Documents/Projects/Sunrise/sps_models/lmcave_ext_noheader.dat', Reader=asciitable.NoHeader)
    extlaminv = extdata['col1']
    ext_AlamAV = extdata['col2']

    #print extlaminv
    #print ext_AlamAV
    extlam_microns = 1.0/extlaminv
    
    #return

    #dust params
    #b = -0.5
    s = 0.35 #0.35 #0.35 #s = 1.35
    b = -0.5

    galscale_kpc = 20.0

    Nfilters = (filter_array.shape)[0]
    Npixels = (1.0*singleslice.shape[0])**2
    i=0
    count=0
    slicecount = 0
    initz = 0.0
    #    for ns in nums:
    #        fn = 'broadband_'+ns+'.fits'
    for fn in bbfiles:
        count = count + 1
        hdulist = pyfits.open(fn)
        #why? -- oh, to skip the addition if the run is still proceeding
        if len(hdulist) < imind+1:
            continue
        
        camdata = ((hdulist)[imind]).data
        filhdu = hdulist[filter_index] ; bbhdu = hdulist[bb_index] ; iqhdu = hdulist[iq_index] ; camparhdu = hdulist[params_index] ; camhdu = hdulist[imind]
        auxhdu = hdulist[aux_index]

        redshift = bbhdu.header.get('redshift')
        if count % mod == 1:
            initz = redshift

        if redshift < startz:
            continue
        #if redshift > endz:
        #    continue

        

        lambda_eff[i,:] = filhdu.data.field('lambda_eff') #; print filhdu.data.field('L_lambda_eff0').shape
        ewidth_lambda[i,:] = filhdu.data.field('ewidth_lambda')
        ewidth_nu[i,:] = filhdu.data.field('ewidth_nu')
        L_lambda_eff[i,:] = filhdu.data.field('L_lambda_eff0')
        AB_mag[i,:] = filhdu.data.field('AB_mag0')
        L_lambda_eff_ns[i,:] = filhdu.data.field('L_lambda_eff_nonscatter0')
        AB_mag_ns[i,:] = filhdu.data.field('AB_mag_nonscatter0')

        z_array[i]=redshift
        pix_array[i]=camhdu.header.get('CD1_1')
        camD_array[i]=camparhdu.header.get('cameradist')
        fov_array[i]=camparhdu.header.get('linear_fov')
        
        L_lambda[i,:]=iqhdu.data.field('L_lambda')
        L_lambda_ns[i,:]=iqhdu.data.field('L_lambda_nonscatter0')
        L_lambda_out[i,:]=iqhdu.data.field('L_lambda_out0')
        L_lambda_abs[i,:]=iqhdu.data.field('L_lambda_absorbed')


        nzindices = np.where(camdata[6,:,:] > 0.0)  #nonzero indices -- this is checking the H-band WFC3 I assume?
        numnz = 1.0*(nzindices[0].shape)[0]
        nzfrac = numnz/Npixels
        nz_kpcsq = nzfrac * pix_array[i]**2


        print fn, redshift #, Npixels**0.5, numnz, nz_kpcsq, nzfrac, np.min(zweightslice), np.max(zweightslice)
        totalflux = totalflux + camdata
        partialflux = partialflux + camdata


        #is this part uber slow??
        zweightslice = np.where(totalflux > 0.0, ( ((totalflux-camdata)*zweightslice + camdata*redshift)/totalflux ), cambasezeros)
        partialzweight = np.where(partialflux > 0.0, ( ((partialflux-camdata)*partialzweight + camdata*redshift)/partialflux ), cambasezeros)

        if has_aux==True:
            auxdata = auxhdu.data
            smdata = (auxdata)[4,:,:] ; smetalsdata = (auxdata)[5,:,:]  # I THINK
            sfrdata = (auxdata)[2,:,:]

            gasmass = (auxdata)[0,:,:] ; metalmass = (auxdata)[1,:,:] ; sfrimage = (auxdata)[2,:,:]  ;  gastemp = (auxdata)[3,:,:]
            restframe_lambda_eff = ( (1.0e6)*lambda_eff[i,:]/(1.0 + redshift) )
            restframe_invlambda_eff = 1.0/restframe_lambda_eff

            Alambda_eff = np.interp(restframe_invlambda_eff.flatten(), extlaminv, ext_AlamAV)

            tau_f = np.zeros_like(gasmass)
            block_dandanxu = np.zeros_like(camdata)
            block_gfs = np.zeros_like(camdata)
            constantblock = np.zeros_like(gasmass)

            nonzero_gasmass_ind = np.where(gasmass > 0.0)[0] #; print nonzero_gasmass_ind.shape
            finite_gasmass_ind = np.isfinite(gasmass)
        
            sigma_pixels = (galscale_kpc/pix_array[i])/(1.0 + redshift) ; print 'sigma in pixels    ', sigma_pixels
            t = scipy.signal.gaussian(101,sigma_pixels)
            kernel = t.reshape(101,1)*t.reshape(1,101)
            kernel /= kernel.sum()

            gasmass_convolved = gasmass #scipy.ndimage.filters.gaussian_filter( gasmass, sigma_pixels)
            metalmass_convolved = metalmass #scipy.ndimage.filters.gaussian_filter( metalmass, sigma_pixels)
            logtemp = np.log10(gastemp)
            covering_factor = 1.0 #0.33 #0.1 + 0.9*(logtemp - 2.5)/(6.1-2.5)   'MOD5'

            metallicity_zsolar = (metalmass_convolved/gasmass_convolved)/0.02
        
            constantblock = covering_factor*((1.0+redshift)**b)*((metallicity_zsolar)**s)*(metalmass_convolved/0.02)*4.28e-8
        
            for f in range(Nfilters):
                if (nonzero_gasmass_ind.shape)[0] < Npixels:
                    continue
                Alam = Alambda_eff[f]
                tau_f = (Alam)*constantblock
                block_dandanxu[f,:,:] = (camdata[f,:,:]/tau_f)*(1.0 - np.exp(-1.0*tau_f))
                block_gfs[f,:,:] = camdata[f,:,:]*np.exp(-1.0*tau_f) 
                if f==6 or f==2:
                    print Alam, np.max(logtemp),np.min(logtemp),np.max(covering_factor),np.min(covering_factor), np.mean(tau_f), np.max(tau_f), np.max(metallicity_zsolar)

        totalflux_dandanxu = totalflux_dandanxu + block_dandanxu
        totalflux_gfs = totalflux_gfs + block_gfs
        #print np.sum(totalflux), np.sum(totalflux_dandanxu)
        intSB[i,:] = np.sum(np.sum(camdata,axis=1),axis=1)
        cumSB[i,:] = np.sum(np.sum(totalflux,axis=1),axis=1)
        

        if has_aux==True:
            totalsm = totalsm + smdata*(pix_array[i]**2)  #want per pixel^2
            partialsm = partialsm + smdata*(pix_array[i]**2)
            totalsfr = totalsfr + sfrdata*(pix_array[i]**2)  #want per pixel^2
            partialsfr = partialsfr + sfrdata*(pix_array[i]**2)
            totalstellarmetals = totalstellarmetals + smetalsdata*(pix_array[i]**2)  #want per pixel^2
            partialstellarmetals = partialstellarmetals + smetalsdata*(pix_array[i]**2)
            sm_array[i] = np.sum(np.sum(smdata*(pix_array[i]**2),axis=0),axis=0)
            mets_array[i] = np.sum(np.sum(smetalsdata*(pix_array[i]**2),axis=0),axis=0)
        
        if (count % mod == 0) or (fn==bbfiles[-1]):
            #write partial flux
            slicecount = slicecount+1
            imagefile = os.path.join(dirname, 'bbslice_{:04d}'.format(slicecount)+'.fits')
            primhdu = pyfits.PrimaryHDU(partialflux)
            primhdu.header.update('zmin',initz)
            primhdu.header.update('zmax',redshift)
            primhdu.update_ext_name('sci') ; primhdu.header.update('IMUNIT',fluxunit)
            #sechdu = pyfits.ImageHDU(totalflux)
            #sechdu.update_ext_name('sci_cumulative')
            sechdu = pyfits.ImageHDU(partialsm)
            sechdu.update_ext_name('STELLAR_MASS') ; sechdu.header.update('IMUNIT',stellarmassunit)
            terhdu = pyfits.ImageHDU(partialstellarmetals)
            terhdu.update_ext_name('STELLAR_METALS_MASS') ; terhdu.header.update('IMUNIT',metalsunit)
            sfrhdu = pyfits.ImageHDU(partialsfr)
            sfrhdu.update_ext_name('SFR') ; sechdu.header.update('IMUNIT',stellarmassunit+'/yr')

            tempzhdu = pyfits.ImageHDU(partialzweight)

            tempzhdu.update_ext_name('ZWEIGHTSLICE')

            newlist = pyfits.HDUList([primhdu,filhdu,sechdu,terhdu,tempzhdu,sfrhdu])

            newlist.writeto(imagefile,clobber=True)
            partialflux = np.zeros_like(cambase) ; partialsm = np.zeros_like(singleslice) ; partialstellarmetals = np.zeros_like(singleslice)
            partialzweight = np.zeros_like(cambase)

        i=i+1
        hdulist.close()



    primhdu = pyfits.PrimaryHDU(totalflux) ; primhdu.header.update('IMUNIT',fluxunit) ; primhdu.header.update('PIXSCALE', pixel_fov_arcsec)
    sfrhdu = pyfits.ImageHDU(totalsfr)
    sfrhdu.update_ext_name('SFR') ; sfrhdu.header.update('IMUNIT',stellarmassunit+'/yr')
    sechdu = pyfits.ImageHDU(totalsm)
    sechdu.update_ext_name('STELLAR_MASS') ; sechdu.header.update('IMUNIT',stellarmassunit)
    terhdu = pyfits.ImageHDU(totalstellarmetals)
    terhdu.update_ext_name('STELLAR_METALS_MASS') ; terhdu.header.update('IMUNIT',metalsunit)
    dusthdu = pyfits.ImageHDU(totalflux_dandanxu) ; dusthdu.header.update('IMUNIT',fluxunit) ; dusthdu.update_ext_name('SEMIANALYTIC_DUST')
    gfshdu = pyfits.ImageHDU(totalflux_gfs) ; gfshdu.header.update('IMUNIT',fluxunit) ; gfshdu.update_ext_name('GFS_SCREEN_DUST')
    
    #zweightvalues = np.where(totalflux > 0.0, zweightslice/totalflux, cambasezeros - 1.0)
    zhdu = pyfits.ImageHDU(zweightslice) ; zhdu.update_ext_name('WEIGHTED_Z')
    
    L_lambda_hdu = pyfits.ImageHDU(L_lambda) ; L_lambda_hdu.update_ext_name('L_lambda') ; L_lambda_hdu.header.update('UNIT',L_lambda_unit)
    L_lambda_ns_hdu = pyfits.ImageHDU(L_lambda_ns) ; L_lambda_ns_hdu.update_ext_name('L_lambda_nonscatter') ; L_lambda_ns_hdu.header.update('UNIT',L_lambda_ns_unit)
    L_lambda_out_hdu = pyfits.ImageHDU(L_lambda_out) ; L_lambda_out_hdu.update_ext_name('L_lambda_out') ; L_lambda_out_hdu.header.update('UNIT',L_lambda_out_unit)
    L_lambda_abs_hdu = pyfits.ImageHDU(L_lambda_abs) ; L_lambda_abs_hdu.update_ext_name('L_lambda_absorption') ; L_lambda_abs_hdu.header.update('UNIT',L_lambda_abs_unit)

    filhdu1 = pyfits.ImageHDU(lambda_eff) ; filhdu1.update_ext_name('lambda_eff') ; filhdu1.header.update('UNIT',lambda_eff_unit)
    filhdu2 = pyfits.ImageHDU(ewidth_lambda) ; filhdu2.update_ext_name('ewidth_lambda') ; filhdu2.header.update('UNIT',ewidth_lambda_unit)
    filhdu3 = pyfits.ImageHDU(ewidth_nu) ; filhdu3.update_ext_name('ewidth_nu') ; filhdu3.header.update('UNIT',ewidth_nu_unit)
    filhdu4 = pyfits.ImageHDU(L_lambda_eff) ; filhdu4.update_ext_name('L_lambda_eff') ; filhdu4.header.update('UNIT',L_lambda_eff_unit)
    filhdu5 = pyfits.ImageHDU(AB_mag) ; filhdu5.update_ext_name('AB_mag')
    filhdu6 = pyfits.ImageHDU(L_lambda_eff_ns) ; filhdu6.update_ext_name('L_lambda_eff_nonscatter') ; filhdu6.header.update('UNIT',L_lambda_eff_ns_unit)
    filhdu7 = pyfits.ImageHDU(AB_mag_ns) ; filhdu7.update_ext_name('AB_mag_nonscatter')
    filhdu8 = pyfits.ImageHDU(intSB) ; filhdu8.update_ext_name('total_SB') ; filhdu8.header.update('UNIT',intSB_unit)
    filhdu9 = pyfits.ImageHDU(cumSB) ; filhdu9.update_ext_name('cumulative_SB') ; filhdu9.header.update('UNIT',cumSB_unit)

    col1 = pyfits.Column(name='redshift',format='D',array=z_array)
    col2 = pyfits.Column(name='pixelscale',format='D',unit='kpc',array=pix_array)
    col3 = pyfits.Column(name='cameradistance',format='D',unit='kpc',array=camD_array)
    col4 = pyfits.Column(name='linearfov',format='D',unit='kpc',array=fov_array)
    col5 = pyfits.Column(name='stellarmass',format='D',unit=stellarmassunit,array=sm_array)
    col6 = pyfits.Column(name='stellarmetalsmass',format='D',unit=metalsunit,array=mets_array)
    cols = pyfits.ColDefs([col1,col2,col3,col4,col5,col6]) ; geometryhdu = pyfits.new_table(cols) ; geometryhdu.update_ext_name('GEOMETRY')

    newlist = pyfits.HDUList([primhdu, dusthdu, gfshdu, sechdu, terhdu, sfrhdu, zhdu, geometryhdu, lamhdu, L_lambda_hdu, L_lambda_ns_hdu, L_lambda_out_hdu, L_lambda_abs_hdu,
                              filname_hdu, filhdu1,filhdu2,filhdu3,filhdu4,filhdu5,filhdu6,filhdu7,filhdu8,filhdu9])
    newlist.writeto(os.path.join(dirname,'bbdata.fits'),clobber=True)


    result = write_filtered_images(totalflux,filter_hdu,dirname=dirname,imunit=fluxunit)


    subprocess.call(['cp','sfrhist_base.stub',dirname])
    subprocess.call(['cp','mcrx_base.stub',dirname])
    subprocess.call(['cp','broadband_base.stub',dirname])
    subprocess.call(['cp','simpar',dirname])

    
    return totalflux



def write_filtered_images(totalflux,filterhdu,dirname='SCI_IMAGES',imunit='W/m/m^2/Sr'):
    print "Writing summed output image for each filter..."
    filterdata = filterhdu.data
    filters = filterdata.field('filter')

        
    for filname in filters:
        splitted = (filname.split('.res'))[0]
        imagefile = os.path.join(dirname,splitted+'.fits')
        index = np.where(filters == filname)
        newimage = totalflux[index[0][0],:,:]
        primhdu = pyfits.PrimaryHDU(newimage)
        primhdu.update_ext_name('sci') ; primhdu.header.update('IMUNIT',imunit)
        #print primhdu.header.ascardlist()
        newlist = pyfits.HDUList([primhdu])
        newlist.writeto(imagefile,clobber=True)
        print imagefile, newimage.shape
    return 1



if __name__=="__main__":
    outputdir = sys.argv[1]
    print "Creating LightCone Images in directory: ", outputdir

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)


    res = add_all_broadband_flux(dirname=outputdir,startz=0.5)

