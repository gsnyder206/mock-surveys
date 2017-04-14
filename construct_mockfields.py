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



#construct real illustris lightcones from individual images

#in parallel, produce estimated Hydro-ART surveys based on matching algorithms -- high-res?


def build_lightcone_images(image_info_file,run_type='images'):

    data=ascii.read(image_info_file)
    print(data)

    full_npix=data['full_npix'][0]

    #get expected shape
    test_file=os.path.join(data['run_dir'][0],'broadbandz.fits')
    tfo =pyfits.open(test_file)
    print(tfo.info())
    cube=tfo['CAMERA0-BROADBAND-NONSCATTER'].data
    cubeshape=cube.shape
    print(cubeshape)
    auxcube=tfo['CAMERA0-AUX'].data

    filters_hdu = tfo['FILTERS']

    lightcone_dir=os.path.abspath(os.path.dirname(image_info_file))
    print('Constructing lightcone data from: ', lightcone_dir)

    output_dir = os.path.join(lightcone_dir,os.path.basename(image_info_file).rstrip('.txt'))
    print('Saving lightcone outputs in: ', output_dir)
    if not os.path.lexists(output_dir):
        os.mkdir(output_dir)

    success_catalog=os.path.join(output_dir,os.path.basename(image_info_file).rstrip('.txt')+'_success.txt')

    image_filelabel='lightcone_image'

    N_filters = cubeshape[0]
    N_aux=auxcube.shape[0]

    image_cube = np.zeros((N_filters,full_npix,full_npix),dtype=np.float64)
    aux_cube = np.zeros((N_aux,full_npix,full_npix),dtype=np.float64)

    success=[]

    #for bigger files, may need to split by filter first

    for origin_i,origin_j,run_dir,this_npix in zip(data['origin_i'],data['origin_j'],data['run_dir'],data['this_npix']):
        try:
            bblist=pyfits.open(os.path.join(run_dir,'broadbandz.fits'))
            this_cube = bblist['CAMERA0-BROADBAND-NONSCATTER'].data
            bblist.close()
            success.append(True)
        except:
            print('Missing file, ', run_dir)
            success.append(False)
            continue
        i_tc=0
        j_tc=0
        i_tc1=this_npix
        j_tc1=this_npix

        if origin_i < 0:
            i0=0
            i_tc=-1*origin_i
        else:
            i0=origin_i

        if origin_j < 0:
            j0=0
            j_tc=-1*origin_j
        else:
            j0=origin_j

        if i0+this_npix > full_npix:
            i1=full_npix
            i_tc1= full_npix-i0   #this_npix - (i0+this_npix-full_npix)
        else:
            i1=i0+this_npix-i_tc

        if j0+this_npix > full_npix:
            j1=full_npix
            j_tc1= full_npix-j0
        else:
            j1=j0+this_npix-j_tc


        sub_cube1=image_cube[:,i0:i1,j0:j1]
        this_subcube=this_cube[:,i_tc:i_tc1,j_tc:j_tc1]
        print(origin_i,origin_j,this_cube.shape, this_npix, sub_cube1.shape, this_subcube.shape)

        image_cube[:,i0:i1,j0:j1] = sub_cube1 + this_subcube


    success=np.asarray(success)
    newcol=astropy.table.column.Column(data=success,name='success')
    data.add_column(newcol)
    ascii.write(data,output=success_catalog)

    filters_data=filters_hdu.data

    for i,filname in enumerate(filters_data['filter']):
        print(filname)
        outname=os.path.join(output_dir,image_filelabel+'_'+filname.replace('/','-')+'.fits')
        print('saving:', outname)

        primary_hdu=pyfits.PrimaryHDU(image_cube[i,:,:])
        output_list=pyfits.HDUList([primary_hdu])
        output_list.writeto(outname,overwrite=True)
        output_list.close()


    #convert units before saving.. or save both?

    return
