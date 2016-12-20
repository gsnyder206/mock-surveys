

####
#### Name:    setup_cylinders.py
#### Author:  Greg Snyder   gsnyder@stsci.edu
#### Purpose:  Converts lightcone parameter file into Sunrise input files (Jonsson 2006).
		i.e., sfrhist*.config, mcrx*.config, and broadband*.config
		to be used in conjunction with *_base.stub files.
#### Disclaimer:  This code is provided AS-IS with absolutely NO warranty.
####           It is largely meant as a guide rather than ideal code.  
####           It can and should be replaced with other lightcone generation techniques. 
####           I make no claims about the immediate usability of this code.
####           That said, I am happy to field questions and discuss issues 
####           related to this code. And also to help use it.
#### License:  ?
#### Credit:   Users should cite Snyder et al. (2014): "Mock UDFs from Hydro Sims with Sunrise"
####           AND ALSO   Kitzbichler & White 2007, Henriques et al. 2013, Overzier et al. 2013, etc.
####           These papers were used as source material to create the default lightcone algorithm.
####  






import math
import string
import sys
import os
import struct
import numpy as np
import asciitable
import pyfits
import subprocess

if __name__=="__main__":
    if (len(sys.argv)!=2):
        print "Error: requires a filename arguments for the lightcone generation master file."
        quit()
    
    masterfile=(sys.argv)[1]

    subprocess.call(['cp',masterfile,'.'])
    
    n='4'      #how many total processors do you want?
    tile='4'   #how many processors do you want to reserve per node?#location="~gsnyder/sunrise/bin_Apr13_2012"
    
    location="/n/home09/gsnyder/sunrise/bin_GFM_May2013_DNDEBUG"  #location of sunrise binary executable files
    queue="sunrise"    #which queue to submit to? try short_parallel for debugging?  Maybe?
    snapexten=".hdf5"      #hack to get the script to recognize different names of simulation input files ".hdf5" for that extension, blank ("") for binary files
    memlimit="12505856"  #in KB   #memlimit=90000000

    snapbase="cylinder_"

    print "Setting up cosmological lightcone Sunrise run... input file (copying here): "+masterfile
    masterdata = asciitable.read(masterfile,comment="#",Reader=asciitable.NoHeader)

    print "WARNING: IT IS THE USERS RESPONSIBILITY TO ENSURE THE COLUMNS NAME MEAN THE RIGHT THING!  If it is much later than Aug 2, 2012, may need to check..."

    f = open(masterfile,'r')
    contents = f.readlines()
    print "Masterfile has length, in lines: ", len(contents)
    u3line = contents[10]
    u1line = contents[11]
    u2line = contents[12]
    print "I think the unit vectors are in lines: "
    print u3line, u1line, u2line
    u3string = ((u3line.split(' [')[1]).split(']')[0]).split() #parse into a list of strings
    u3vec = map(float,np.array(u3string))
    print "I parsed the unit vector to be: ", u3vec

    u2string = ((u2line.split(' [')[1]).split(']')[0]).split() #parse into a list of strings
    u2vec = map(float,np.array(u2string))
    print "I parsed the up vector to be: ", u2vec

    id_int = masterdata['col1']  #I think these all end up being strings
    v_In_phys_x = masterdata['col10']
    v_In_phys_y = masterdata['col11']
    v_In_phys_z = masterdata['col12']

    cam_phys_x = masterdata['col13']
    cam_phys_y = masterdata['col14']
    cam_phys_z = masterdata['col15']

    cam_phys_offset_x = masterdata['col16']
    cam_phys_offset_y = masterdata['col17']
    cam_phys_offset_z = masterdata['col18']

    cam_distance_phys = (cam_phys_offset_x**2 + cam_phys_offset_y**2 + cam_phys_offset_z**2)**0.5

    fov_kpc = masterdata['col19']
    box_redshift = masterdata['col20']

    fov_rad = 2.0*np.arctan(fov_kpc/(2.0*cam_distance_phys))
    amBline = contents[8]
    amstring = (amBline.split(':'))[1]
    fov_arcmin_file = float(amstring)
    fov_rad_file = (math.pi/180.0)*(fov_arcmin_file/60.0)

    print "I parsed the field of view, in radians, along the small direction, to be: ",fov_rad_file
    
    Ncylinders = int((fov_kpc.shape)[0])
    print '  Setting up for '+str(Ncylinders)+' cylinders...'

    for id in id_int[1:]:  #really want to skip that first one...
        index = (np.where(id_int==id))[0]
        numstring = '{:04d}'.format(id)
        snapfile=snapbase+numstring+snapexten
        exists = os.path.exists(snapfile)
        if exists==False:
            print snapfile+' does not exist, skipping...'
            continue
        #assert exists, "Ack! "+snapfile+" does not exist."
        bsubfile="sunrise_"+numstring+'.bsub'
        sfile='sfrhist_'+numstring+'.config'
        mfile='mcrx_'+numstring+'.config'
        bfile='broadband_'+numstring+'.config'

        bsubf = open(bsubfile,'w')
        bsubf.write('#!/bin/bash\n')
        bsubf.write('#BSUB -q '+queue+'\n')
        bsubf.write('#BSUB -J cylinder_'+numstring+'\n')
        bsubf.write('#BSUB -n '+n+'\n')
        bsubf.write('#BSUB -oo cylinder_'+numstring+'_lsf.out\n')
        bsubf.write('#BSUB -eo cylinder_'+numstring+'_lsf.err\n')
        bsubf.write('#BSUB -M '+memlimit+'\n')
        bsubf.write('#BSUB -R \"span[ptile='+tile+']\"\n')
        bsubf.write('#BSUB -g /GFS2\n')
        bsubf.write('#BSUB -C 0\n\n')
        bsubf.write(location+'/sfrhist '+sfile+' 1> sfrhist_'+numstring+'.out 2>&1\n')
        bsubf.write(location+'/mcrx '+mfile+' 1> mcrx_'+numstring+'.out 2>&1\n')
        bsubf.write(location+'/broadband '+bfile+' 1> broadband_'+numstring+'.out 2>&1\n')
        bsubf.write('rm -f mcrx_'+numstring+'.fits\n')
        bsubf.write('rm -f sfrhist_'+numstring+'.fits\n')
        bsubf.close()


        sf = open(sfile,'w')
        sf.write('#Lightcone Cylinder Parameter File for Sunrise, sfrhist\n')
        sf.write('include_file\t\t sfrhist_base.stub\n')
        sf.write('snapshot_file\t\t '+snapfile+'\n')
        sf.write('output_file\t\t sfrhist_'+numstring+'.fits\n')
        sf.write('translate_origin\t {:15.3f} {:15.3f} {:15.3f}\n'.format(float(v_In_phys_x[index]),float(v_In_phys_y[index]),float(v_In_phys_z[index]))+'\n')
        sf.close()

        mf = open(mfile,'w')
        mf.write('#Lightcone Cylinder Parameter File for Sunrise, mcrx\n')
        mf.write('include_file\t\t mcrx_base.stub\n')
        mf.write('input_file\t\t sfrhist_'+numstring+'.fits\n')
        mf.write('output_file\t\t mcrx_'+numstring+'.fits\n')
        mf.write('camera_position\t\t {:15.3f} {:15.3f} {:15.3f}\n'.format(float(cam_phys_offset_x[index]),float(cam_phys_offset_y[index]),float(cam_phys_offset_z[index])))
        mf.write('camera_direction\t {:15.8f} {:15.8f} {:15.8f}\n'.format(u3vec[0],u3vec[1],u3vec[2]))
        mf.write('camera_up\t\t {:15.8f} {:15.8f} {:15.8f}\n'.format(u2vec[0],u2vec[1],u2vec[2]))
        mf.write('camerafov\t\t {:15.10f}\n'.format(float(fov_rad[index])))
        #mf.write('translate_origin\t {:15.3f} {:15.3f} {:15.3f}\n'.format(float(v_In_phys_x[index]),float(v_In_phys_y[index]),float(v_In_phys_z[index]))+'\n')
        mf.close()

        bf = open(bfile,'w')
        bf.write('#Lightcone Cylinder Parameter File for Sunrise, broadband\n')
        bf.write('include_file\t\t broadband_base.stub\n')
        bf.write('input_file\t\t mcrx_'+numstring+'.fits\n')
        bf.write('output_file\t\t broadband_'+numstring+'.fits\n')
        bf.write('redshift\t\t {:10.6f}\n'.format(float(box_redshift[index])))
        bf.close()


        print snapfile,bsubfile,sfile,mfile,bfile,box_redshift[index],fov_rad[index]
    





