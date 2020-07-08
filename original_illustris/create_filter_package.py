import numpy as np
import astropy
import astropy.cosmology
import astropy.io.fits as pyfits
import astropy.units as u
import astropy.io.ascii as ascii

import os
import pandas



def hdu_from_existing_filter(input_file):
    return hdu


def wfirst_filters():
    
    xls='WFIRST_WIMWSM_throughput_data_190531.xlsm'
    
    throughputs=pandas.read_excel(xls,'EffectiveArea')

    wl_microns=throughputs['Unnamed: 0'][19:].values

    r062_A=throughputs['Unnamed: 1'][19:].values
    z087_A=throughputs['Unnamed: 2'][19:].values
    y106_A=throughputs['Unnamed: 3'][19:].values
    j129_A=throughputs['Unnamed: 4'][19:].values
    w146_A=throughputs['Unnamed: 5'][19:].values
    h158_A=throughputs['Unnamed: 6'][19:].values
    f184_A=throughputs['Unnamed: 7'][19:].values

    maxA = np.max(throughputs[['Unnamed: 2','Unnamed: 3','Unnamed: 4',
                               'Unnamed: 5','Unnamed: 6','Unnamed: 7','Unnamed: 1']][19:].values)
    
    r_tr = r062_A/maxA
    z_tr = z087_A/maxA
    y_tr = y106_A/maxA
    j_tr = j129_A/maxA
    w_tr = w146_A/maxA
    h_tr = h158_A/maxA
    f_tr = f184_A/maxA
    
    rhdu = wf_hdu(r_tr[r_tr > 1.0e-4],wl_microns[r_tr > 1.0e-4])
    rhdu.header['EXTNAME']='wfirst/wfi_r062'
    
    zhdu = wf_hdu(z_tr[z_tr > 1.0e-4],wl_microns[z_tr > 1.0e-4])
    zhdu.header['EXTNAME']='wfirst/wfi_z087'
    
    yhdu = wf_hdu(y_tr[y_tr > 1.0e-4],wl_microns[y_tr > 1.0e-4])
    yhdu.header['EXTNAME']='wfirst/wfi_y106'
    
    jhdu = wf_hdu(j_tr[j_tr > 1.0e-4],wl_microns[j_tr > 1.0e-4])
    jhdu.header['EXTNAME']='wfirst/wfi_j129'
    
    whdu = wf_hdu(w_tr[w_tr > 1.0e-4],wl_microns[w_tr > 1.0e-4])
    whdu.header['EXTNAME']='wfirst/wfi_w146'
    
    hhdu = wf_hdu(h_tr[h_tr > 1.0e-4],wl_microns[h_tr > 1.0e-4])
    hhdu.header['EXTNAME']='wfirst/wfi_h158'
    
    fhdu = wf_hdu(f_tr[f_tr > 1.0e-4],wl_microns[f_tr > 1.0e-4])
    fhdu.header['EXTNAME']='wfirst/wfi_f184'
    
    '''
    Z_tr = Z_Aeff_m2/np.max(Z_Aeff_m2)
    Y_tr = Y_Aeff_m2/np.max(Y_Aeff_m2)
    J_tr = J_Aeff_m2/np.max(J_Aeff_m2)
    W_tr = W_Aeff_m2/np.max(W_Aeff_m2)
    H_tr = H_Aeff_m2/np.max(H_Aeff_m2)
    F_tr = F_Aeff_m2/np.max(F_Aeff_m2)

    zhdu = wf_hdu(Z_tr,Z_wl_um)
    zhdu.header['EXTNAME']='WFI_DRM15_Z087'
    yhdu = wf_hdu(Y_tr,Y_wl_um)
    yhdu.header['EXTNAME']='WFI_DRM15_Y106'
    jhdu = wf_hdu(J_tr,J_wl_um)
    jhdu.header['EXTNAME']='WFI_DRM15_J129'
    whdu = wf_hdu(W_tr,W_wl_um)
    whdu.header['EXTNAME']='WFI_DRM15_W149'
    hhdu = wf_hdu(H_tr,H_wl_um)
    hhdu.header['EXTNAME']='WFI_DRM15_H158'
    fhdu = wf_hdu(F_tr,F_wl_um)
    fhdu.header['EXTNAME']='WFI_DRM15_F184'
    '''


    return (rhdu,zhdu,yhdu,jhdu,whdu,hhdu,fhdu)


    
def wfirst_filters_drm15():
    Z_wl_um = np.asarray([0.75,0.7559,0.7617,0.7676,0.7734,0.7793,0.7851,0.7910,0.7968,0.8027,0.8085,0.8144,0.8202,0.8261,0.8320,0.8378,0.8437,0.8495,0.8554,0.8612,0.8671,0.8729,0.8788,0.8846,0.8905,0.8963,0.9022,0.9080,0.9139,0.9198,0.9256,0.9315,0.9373,0.9432,0.9490,0.9549,0.9607,0.9666,0.9724,0.9783,0.9841])
    Z_Aeff_m2 = np.asarray([0.069,1.955,2.317,2.292,2.292,2.294,2.279,2.275,2.266,2.271,2.268,2.272,2.266,2.270,2.262,2.258,2.263,2.247,2.251,2.252,2.273,2.270,2.267,2.265,2.264,2.253,2.258,2.255,2.255,2.253,2.253,2.249,2.246,2.235,2.231,2.219,2.214,2.117,0.050,0.000,0.000])
    Y_wl_um = np.asarray([0.92,0.9280,0.9361,0.9441,0.9522,0.9602,0.9683,0.9763,0.9844,0.9924,1.0005,1.0085,1.0166,1.0246,1.0327,1.0407,1.0488,1.0568,1.0649,1.0729,1.0810,1.0890,1.0971,1.1051,1.1132,1.1212,1.1293,1.1373,1.1454,1.1534,1.1615,1.1695,1.1776,1.1856,1.1937,1.2017,1.2098,1.2178,1.2259,1.2339,1.2420,1.2500])
    Y_Aeff_m2 = np.asarray([0.204,0.743,1.282,1.819,2.139,2.138,2.137,2.136,2.135,2.134,2.133,2.138,2.147,2.153,2.162,2.168,2.177,2.184,2.194,2.205,2.215,2.225,2.233,2.244,2.254,2.264,2.273,2.286,2.294,2.306,2.319,2.328,2.242,1.791,1.338,0.879,0.883,0.418,0.000,0.000,0.000,0.000])
    J_wl_um = np.asarray([1.1,1.1098,1.1195,1.1293,1.1390,1.1488,1.1585,1.1683,1.1780,1.1878,1.1976,1.2073,1.2171,1.2266,1.2366,1.2463,1.2561,1.2659,1.2756,1.2854,1.2951,1.3049,1.3146,1.3244,1.3341,1.3439,1.3537,1.3634,1.3732,1.3829,1.3927,1.4024,1.4122,1.4220,1.4317,1.4415,1.4512,1.4610,1.4707,1.4805,1.4902,1.5000])
    J_Aeff_m2 = np.asarray([0.000,0.000,0.208,0.678,1.154,1.637,2.124,2.328,2.341,2.354,2.368,2.380,2.393,2.406,2.420,2.433,2.446,2.462,2.475,2.487,2.500,2.512,2.524,2.536,2.548,2.559,2.570,2.581,2.592,2.607,2.618,2.626,2.633,2.643,2.430,2.012,1.590,1.164,0.736,0.304,0.000,0.000])
    W_wl_um = np.asarray([0.8,0.8293,0.8585,0.8878,0.9171,0.9463,0.9756,1.0049,1.0341,1.0634,1.0927,1.1220,1.1512,1.1805,1.2098,1.2390,1.2683,1.2976,1.3268,1.3561,1.3854,1.4146,1.4439,1.4732,1.5024,1.5317,1.5610,1.5902,1.6195,1.6488,1.6780,1.7073,1.7366,1.7659,1.7951,1.8247,1.8537,1.8829,1.9122,1.9412,1.9707,2.0000])
    W_Aeff_m2 = np.asarray([0.000,0.000,0.000,0.000,0.153,1.864,1.871,1.948,1.987,2.052,2.087,1.982,2.147,2.226,2.263,2.221,2.299,2.345,2.346,2.429,2.270,2.358,2.444,2.460,2.510,2.530,2.555,2.579,2.576,2.650,2.603,2.653,2.639,2.683,2.640,2.648,2.641,2.586,2.651,2.635,2.618,0.002])
    H_wl_um = np.asarray([1.35,1.3622,1.3744,1.3866,1.3988,1.4110,1.4232,1.4354,1.4476,1.4598,1.4720,1.4841,1.4963,1.5085,1.5207,1.5329,1.5451,1.5573,1.5695,1.5817,1.5939,1.6061,1.6183,1.6305,1.6427,1.6549,1.6671,1.6793,1.6915,1.7037,1.7159,1.7280,1.7402,1.7524,1.7646,1.7768,1.7890,1.8012,1.8134,1.8256,1.8378,1.8500])
    H_Aeff_m2 = np.asarray([0.000,0.495,0.937,1.383,1.834,2.630,2.644,2.661,2.674,2.687,2.698,2.707,2.717,2.725,2.735,2.745,2.754,2.762,2.768,2.776,2.784,2.792,2.799,2.806,2.811,2.816,2.823,2.829,2.834,2.840,2.845,2.850,2.779,2.408,2.035,1.660,1.285,0.530,0.152,0.000,0.000,0.000])
    F_wl_um = np.asarray([1.65,1.6610,1.6720,1.6829,1.6939,1.7049,1.7159,1.7268,1.7378,1.7488,1.7598,1.7707,1.7817,1.7927,1.8037,1.8146,1.8256,1.8366,1.8476,1.8585,1.8695,1.8805,1.8915,1.9024,1.9134,1.9244,1.9354,1.9463,1.9573,1.9683,1.9793,1.9902,2.0012,2.0122,2.0232,2.0341,2.0451,2.0561,2.0671,2.0780,2.0890,2.1000])
    F_Aeff_m2 = np.asarray([0.000,0.521,0.877,1.234,1.592,1.952,2.312,2.580,2.585,2.588,2.592,2.595,2.598,2.601,2.601,2.601,2.603,2.605,2.606,2.609,2.611,2.614,2.615,2.616,2.618,2.621,2.623,2.630,2.635,2.618,2.313,1.697,1.390,1.083,0.775,0.466,0.158,0.000,0.000,0.000,0.000,0.000])


    Z_tr = Z_Aeff_m2/np.max(Z_Aeff_m2)
    Y_tr = Y_Aeff_m2/np.max(Y_Aeff_m2)
    J_tr = J_Aeff_m2/np.max(J_Aeff_m2)
    W_tr = W_Aeff_m2/np.max(W_Aeff_m2)
    H_tr = H_Aeff_m2/np.max(H_Aeff_m2)
    F_tr = F_Aeff_m2/np.max(F_Aeff_m2)

    zhdu = wf_hdu(Z_tr,Z_wl_um)
    zhdu.header['EXTNAME']='WFI_DRM15_Z087'
    yhdu = wf_hdu(Y_tr,Y_wl_um)
    yhdu.header['EXTNAME']='WFI_DRM15_Y106'
    jhdu = wf_hdu(J_tr,J_wl_um)
    jhdu.header['EXTNAME']='WFI_DRM15_J129'
    whdu = wf_hdu(W_tr,W_wl_um)
    whdu.header['EXTNAME']='WFI_DRM15_W149'
    hhdu = wf_hdu(H_tr,H_wl_um)
    hhdu.header['EXTNAME']='WFI_DRM15_H158'
    fhdu = wf_hdu(F_tr,F_wl_um)
    fhdu.header['EXTNAME']='WFI_DRM15_F184'



    return (zhdu,yhdu,jhdu,whdu,hhdu,fhdu)

def wf_hdu(tr,um):
    print(np.asarray(tr))
    print(um)
    
    col2 = pyfits.Column(name='totalrelativeefficiency',format='D',unit='relative fraction',array=np.asarray(tr,dtype=np.float64))
    col1 = pyfits.Column(name='wavelength',format='D',unit='angstrom',array=np.asarray(um,dtype=np.float64)*1.0e4)
    cols = pyfits.ColDefs([col1,col2])
    zhdu = pyfits.BinTableHDU.from_columns(cols)
    return zhdu


def output_sunrise_filter_directory(fitsfile,dirname):
    f = pyfits.open(fitsfile)

    if not os.path.lexists(dirname):
        os.mkdir(dirname)

    all_list = []
    grism_list = []
    for hdu in f[1:]:
        #print hdu.header['EXTNAME']
        filter_name = os.path.join(dirname,hdu.header['EXTNAME'])
        filter_wl = hdu.data['wavelength']
        filter_tr = hdu.data['totalrelativeefficiency']
        data = astropy.table.Table([np.round(filter_wl,4),np.round(filter_tr,4)],names=['wl','tr'])
        dn = os.path.dirname(filter_name)
        
        if not os.path.lexists(dn):
            os.mkdir(dn)

        all_list.append(hdu.header['EXTNAME'])

        if hdu.header['EXTNAME'][0:5]=='grism':
            grism_list.append(hdu.header['EXTNAME'])

        ascii.write(data,filter_name,Writer=ascii.NoHeader)



    #output_lists
    #print all_list
    #print grism_list

    st_list = ['hst/wfc3_f336w',
               'hst/acs_f435w',
               'hst/acs_f606w',
               'hst/acs_f814w',
               'hst/wfc3_f125w',
               'hst/wfc3_f160w',
               'wfirst/wfi_r062',
               'wfirst/wfi_z087',
               'wfirst/wfi_y106',
               'wfirst/wfi_j129',
               'wfirst/wfi_w146',
               'wfirst/wfi_h158',
               'wfirst/wfi_f184',
               'jwst/nircam_f115w',
               'jwst/nircam_f150w',
               'jwst/nircam_f200w',
               'jwst/nircam_f277w',
               'jwst/nircam_f356w',
               'jwst/nircam_f444w',
               'jwst/miri_F770W',
               'jwst/miri_F1500W']

    rest_list = ['standard/wfcam_j',
                 'standard/wfcam_h',
                 'standard/wfcam_k',
                 'standard/bessel_b',
                 'standard/bessel_i',
                 'standard/bessel_r',
                 'standard/bessel_u',
                 'standard/bessel_v']

    other_list = ['galex/fuv_1500',
                  'galex/nuv_2500',
                  'sdss/u',
                  'sdss/g',
                  'sdss/r',
                  'sdss/i',
                  'sdss/z',
                  'newfirm/newfirm_h1ccd',
                  'newfirm/newfirm_h2ccd',
                  'newfirm/newfirm_j1ccd',
                  'newfirm/newfirm_j2ccd',
                  'newfirm/newfirm_j3ccd',
                  'newfirm/newfirm_Kccd',
                  'lsst/lsst_u',
                  'lsst/lsst_g',
                  'lsst/lsst_r',
                  'lsst/lsst_i',
                  'lsst/lsst_z',
                  'lsst/lsst_y3']

    all_table = astropy.table.Table([np.asarray(all_list)],names=['filtername'])
    grism_table = astropy.table.Table([np.asarray(grism_list)],names=['filtername'])

    st_table = astropy.table.Table([np.asarray(st_list)],names=['filtername'])
    rest_table = astropy.table.Table([np.asarray(rest_list)],names=['filtername'])
    other_table = astropy.table.Table([np.asarray(other_list)],names=['filtername'])


    ascii.write(all_table,os.path.join(dirname,'filters_all'),Writer=ascii.NoHeader)
    ascii.write(grism_table,os.path.join(dirname,'filters_grism'),Writer=ascii.NoHeader)
    ascii.write(st_table,os.path.join(dirname,'filters_st'),Writer=ascii.NoHeader)
    ascii.write(rest_table,os.path.join(dirname,'filters_rest'),Writer=ascii.NoHeader)
    ascii.write(other_table,os.path.join(dirname,'filters_other'),Writer=ascii.NoHeader)
    

    return


def input_old_filters(flist,folder):
    names = np.asarray( ascii.read(flist,Reader=ascii.NoHeader))
    hdus = []
    #print flist
    #print folder
    #print names


    for n in names:
        #print n[0]
        path = os.path.join(folder,n[0])

        data = ascii.read(path)
        wl_ang = np.round(np.asarray(data['col1']),4)
        tr = np.round(np.asarray(data['col2']),8)
        hdu = wf_hdu(tr,wl_ang/1.0e4)
        hdu.header['EXTNAME']=n[0][0:-4]
        #print hdu.header['EXTNAME']
        hdus.append(hdu)

    return hdus


def miri_filters():
    flist = ['F560W_throughput.fits','F770W_throughput.fits','F1000W_throughput.fits','F1130W_throughput.fits','F1280W_throughput.fits','F1500W_throughput.fits','F1800W_throughput.fits','F2100W_throughput.fits','F2550W_throughput.fits']
    hdus = []

    for file in flist:
        f = pyfits.open(os.path.join('MIRI/filters',file))
        tr = f[1].data['THROUGHPUT']
        wl = f[1].data['WAVELENGTH']
        hdu = wf_hdu(tr,wl/1.0e4)
        hdu.header['EXTNAME']='jwst/miri_'+file[:-16]
        hdus.append(hdu)
    
    return hdus



def nircam_filters():
    flist = ['F070W_throughput.fits',
             'F090W_throughput.fits',
             'F115W_throughput.fits',
             'F150W_throughput.fits',
             'F200W_throughput.fits',
             'F277W_throughput.fits',
             'F356W_throughput.fits',
             'F444W_throughput.fits',
             'F410M_throughput.fits']
    
    hdus = []

    for file in flist:
        f = pyfits.open(os.path.join('NIRCam_Filters',file))
        tr = f[1].data['THROUGHPUT']
        wl = f[1].data['WAVELENGTH']
        hdu = wf_hdu(tr,wl/1.0e4)
        hdu.header['EXTNAME']='jwst/nircam_'+file[:-16]
        hdus.append(hdu)
    
    return hdus


def make_grism_filters(start,stop,number):

    micron_bins = np.logspace(np.log10(start),np.log10(stop),number)
    hdus = []

    for i,b in enumerate(micron_bins[1:]):
        bin_start=micron_bins[i]
        bin_stop=b
        bin_width = bin_stop-bin_start
        bin_center = (bin_start + bin_stop)/2.0
        thisname='grism/gen1_'+'{:5.3f}'.format(bin_center)
        #print bin_start, bin_stop, bin_width, bin_center, thisname
        wl_start = bin_start-bin_width
        wl_stop = bin_stop+bin_width
        wl_grid = np.linspace(wl_start,wl_stop,30.0)
        wl_tr = np.where(np.logical_and(wl_grid >= bin_start,wl_grid < bin_stop ),np.ones_like(wl_grid),np.zeros_like(wl_grid) )
        this_hdu = wf_hdu(wl_tr,wl_grid)
        this_hdu.header['EXTNAME']=thisname
        hdus.append(this_hdu)


    return hdus

if __name__=="__main__":
    #get existing filters and add new ones
    #note these filters are for simple BROADBAND test images, not full simulations, so accuracy=== MEH
    wfirst_hdus = wfirst_filters()
    #print wfirst_hdus[4].data['totalrelativeefficiency'], wfirst_hdus[4].data['wavelength'], wfirst_hdus[4].header.cards

    primhdu = pyfits.PrimaryHDU()
    newlist = pyfits.HDUList([primhdu])
    for hdu in wfirst_hdus:
        newlist.append(hdu)


    all_filters = input_old_filters(os.path.expandvars('$HOME/Dropbox/Projects/FILTERS/sunrise_data/filters_redshifted'),os.path.expandvars('$HOME/Dropbox/Projects/FILTERS/sunrise_data/filters'))
    rest_filters = input_old_filters(os.path.expandvars('$HOME/Dropbox/Projects/FILTERS/sunrise_data/filters_restframe_new'),os.path.expandvars('$HOME/Dropbox/Projects/FILTERS/sunrise_data/filters'))


    for hdu in all_filters:
        newlist.append(hdu)

    for hdu in rest_filters:
        newlist.append(hdu)



    miri_hdus = miri_filters()
    for hdu in miri_hdus:
        newlist.append(hdu)

    nircam_hdus = nircam_filters()
    for hdu in nircam_hdus:
        newlist.append(hdu)

    gen1_grisms = make_grism_filters(0.69,5.0,250)
    #grs_grisms = make_grism_filters(1.35,1.89,250)
    for hdu in gen1_grisms:
        newlist.append(hdu)



    packagef = 'filterpackage.fits'
    newlist.writeto(packagef,clobber=True)

    
    output_sunrise_filter_directory(packagef,'sunrise_filters')
