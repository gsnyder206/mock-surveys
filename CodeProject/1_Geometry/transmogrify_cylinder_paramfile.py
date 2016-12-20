

####
#### Name:    transmogrify_cylinder_paramfile.py
#### Author:  Greg Snyder   gsnyder@stsci.edu
#### Purpose:  Alters parameter setup file for mock survey fields from 
####           continuous-volume hydrodynamical simulations using Sunrise  (Jonsson 2006).
####           The alteration shifts the comoving space in such a way as to view
####           a region of space at different observation epochs.
#### Disclaimer:  This code is provided AS-IS with absolutely NO warranty.
####           It is largely meant as a guide rather than ideal code.  
####           It can and should be replaced with other lightcone generation techniques. 
####           I make no claims about the immediate usability of this code.
####           That said, I am happy to field questions and discuss issues 
####           related to this code. And also to help use it.
#### License:  ?
#### Credit:   Users should cite Snyder et al. (2014): "Mock UDFs from Hydro Sims with Sunrise"
####           AND ALSO   Kitzbichler & White 2007, Henriques et al. 2013, Overzier et al. 2013, etc.
####           These papers were used as source material to create the lightcone generation algorithm.
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



def transmogrify(masterfile,index,outfile):
    print "Transmogrifying lightcone parameter file: "+masterfile
    print "Referencing new lightcone to index: ", index

    masterdata = asciitable.read(masterfile,comment="#",Reader=asciitable.NoHeader)

    print "WARNING: IT IS THE USERS RESPONSIBILITY TO ENSURE THE COLUMNS NAME MEAN THE RIGHT THING!  If it is much later than Aug 2, 2012, may need to check..."

    f = open(masterfile,'r')
    contents = f.readlines()
    print "Masterfile has length, in lines: ", len(contents)

    headerlist = contents[0:35]
    print headerlist
    print "Header has length, in lines: ", len(headerlist)

    cylinderlist = contents[35:]

    
    #u3line = contents[10]
    #u1line = contents[11]
    #u2line = contents[12]
    #print "I think the unit vectors are in lines: "
    #print u3line, u1line, u2line
    #u3string = ((u3line.split(' [')[1]).split(']')[0]).split() #parse into a list of strings
    #u3vec = map(float,np.array(u3string))
    #print "I parsed the unit vector to be: ", u3vec

    #u2string = ((u2line.split(' [')[1]).split(']')[0]).split() #parse into a list of strings
    #u2vec = map(float,np.array(u2string))
    #print "I parsed the up vector to be: ", u2vec

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

    #cam_distance_phys = (cam_phys_offset_x**2 + cam_phys_offset_y**2 + cam_phys_offset_z**2)**0.5

    fov_kpc = masterdata['col19']
    box_redshift = masterdata['col20']
    cylrad = masterdata['col21']


    snaplabel = masterdata['col2']
    snapz = masterdata['col3']
    v_In_co_x = masterdata['col4']
    v_In_co_y = masterdata['col5']
    v_In_co_z = masterdata['col6']
    v_Eg_co_x = masterdata['col7']
    v_Eg_co_y = masterdata['col8']
    v_Eg_co_z = masterdata['col9']

    #fov_rad = 2.0*np.arctan(fov_kpc/(2.0*cam_distance_phys))
    #amBline = contents[8]
    #amstring = (amBline.split(':'))[1]
    #fov_arcmin_file = float(amstring)
    #fov_rad_file = (math.pi/180.0)*(fov_arcmin_file/60.0)

    #print "I parsed the field of view, in radians, along the small direction, to be: ",fov_rad_file
    
    Ncylinders = int((fov_kpc.shape)[0])
    print 'Original paramfile had '+str(Ncylinders)+' cylinders.  Should match:  ', len(cylinderlist)

    print "Referencing new lightcone to index: ", index
    print "Original line entry: "
    print cylinderlist[index]
    ZREF = box_redshift[index]
    print "Reference Redshift: ", ZREF
    AREF = 1.0/(1.0 + ZREF)




    TRANS_id = np.zeros_like(id_int) #shifted
    TRANS_snaplabel = snaplabel  #fixed
    TRANS_snapz = snapz #fixed
    TRANS_v_In_co_x = np.zeros_like(v_In_co_x) #shifted
    TRANS_v_In_co_y = np.zeros_like(v_In_co_x) #shifted
    TRANS_v_In_co_z = np.zeros_like(v_In_co_x) #shifted
    TRANS_v_Eg_co_x = np.zeros_like(v_In_co_x) #shifted
    TRANS_v_Eg_co_y = np.zeros_like(v_In_co_x) #shifted
    TRANS_v_Eg_co_z = np.zeros_like(v_In_co_x) #shifted
    TRANS_v_In_phys_x = np.zeros_like(v_In_co_x) #shifted and modded
    TRANS_v_In_phys_y = np.zeros_like(v_In_co_x) #shifted and modded
    TRANS_v_In_phys_z = np.zeros_like(v_In_co_x) #shifted and modded
    TRANS_cam_phys_x = np.zeros_like(v_In_co_x) #shifted and modded
    TRANS_cam_phys_y = np.zeros_like(v_In_co_x) #shifted and modded
    TRANS_cam_phys_z = np.zeros_like(v_In_co_x) #shifted and modded
    TRANS_cam_phys_offset_x = np.zeros_like(v_In_co_x) #shifted and modded
    TRANS_cam_phys_offset_y = np.zeros_like(v_In_co_x) #shifted and modded
    TRANS_cam_phys_offset_z = np.zeros_like(v_In_co_x) #shifted and modded
    TRANS_fov_kpc = np.zeros_like(v_In_co_x)  #shifted and modded
    TRANS_box_redshift = np.zeros_like(v_In_co_x) #modded
    TRANS_cylrad = np.zeros_like(cylrad)  #shifted


    for i in range(Ncylinders):

        ORIGZCEN = box_redshift[i]

        #the expansion factor at the cosmic time of this section of the lightcone
        AORIG = 1.0/(1.0 + ORIGZCEN)

        #new relative redshift, from a1/a2 = 1 + z'
        NEWZCEN = (1.0 + ORIGZCEN)/(1.0 + ZREF) - 1.0



        ishifted=0
        if i >= index:
            ishifted = i - index

        #to put physical quantities into comoving lightcone quantities
        ZSHIFTED = box_redshift[ishifted]
        ASHIFTED = 1.0/(1.0 + ZSHIFTED)
        AFACTOR = AORIG/ASHIFTED
        #print ishifted, AFACTOR, (ishifted), (AFACTOR).shape, (cam_phys_y[ishifted]).shape
        
        #print ishifted, TRANS_snaplabel[i], TRANS_snapz[i], v_In_co_x[ishifted], v_In_co_y[ishifted], v_In_co_z[ishifted], v_Eg_co_x[ishifted], v_Eg_co_y[ishifted], v_Eg_co_z[ishifted], AFACTOR*v_In_phys_x[ishifted], AFACTOR*v_In_phys_y[ishifted], AFACTOR*v_In_phys_z[ishifted], AFACTOR*cam_phys_x[ishifted], AFACTOR*cam_phys_y[ishifted], AFACTOR*cam_phys_z[ishifted], AFACTOR*cam_phys_offset_x[ishifted], AFACTOR*cam_phys_offset_y[ishifted], AFACTOR*cam_phys_offset_z[ishifted], AFACTOR*fov_kpc[ishifted], NEWZCEN, cylrad[ishifted]

        TRANS_id[i] = ishifted
        TRANS_v_In_co_x[i] = v_In_co_x[ishifted] #shifted
        TRANS_v_In_co_y[i] = v_In_co_y[ishifted] #shifted
        TRANS_v_In_co_z[i] = v_In_co_z[ishifted] #shifted
        TRANS_v_Eg_co_x[i] = v_Eg_co_x[ishifted] #shifted
        TRANS_v_Eg_co_y[i] = v_Eg_co_y[ishifted] #shifted
        TRANS_v_Eg_co_z[i] = v_Eg_co_z[ishifted] #shifted
        TRANS_v_In_phys_x[i] = AFACTOR*v_In_phys_x[ishifted] #shifted and modded
        TRANS_v_In_phys_y[i] = AFACTOR*v_In_phys_y[ishifted] #shifted and modded
        TRANS_v_In_phys_z[i] = AFACTOR*v_In_phys_z[ishifted] #shifted and modded
        TRANS_cam_phys_x[i] =  AFACTOR*cam_phys_x[ishifted] #shifted and modded
        TRANS_cam_phys_y[i] = AFACTOR*cam_phys_y[ishifted] #shifted and modded
        TRANS_cam_phys_z[i] = AFACTOR*cam_phys_z[ishifted] #shifted and modded
        TRANS_cam_phys_offset_x[i] = AFACTOR*cam_phys_offset_x[ishifted] #shifted and modded
        TRANS_cam_phys_offset_y[i] = AFACTOR*cam_phys_offset_y[ishifted] #shifted and modded
        TRANS_cam_phys_offset_z[i] =  AFACTOR*cam_phys_offset_z[ishifted] #shifted and modded
        TRANS_fov_kpc[i] =  AFACTOR*fov_kpc[ishifted] #shifted and modded
        TRANS_box_redshift[i] = NEWZCEN #modded
        TRANS_cylrad[i] = cylrad[ishifted] #shifted






    mogfile = open(outfile,'w')
    for hline in headerlist:
        mogfile.write(hline) #; print hline



    Nnew = Ncylinders - index
    print Nnew

    for i in np.arange(Nnew)+index:
        line = "{:5d}   {:4s}   {:7.4f}   {:10.4f}   {:10.4f}   {:10.4f}" \
            "   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}" \
            "   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}" \
            "   {:10.4f}   {:7.4f}   {:10.4f}\n".format(TRANS_id[i],TRANS_snaplabel[i],TRANS_snapz[i],
            TRANS_v_In_co_x[i], TRANS_v_In_co_y[i], TRANS_v_In_co_z[i], 
            TRANS_v_Eg_co_x[i], TRANS_v_Eg_co_y[i], TRANS_v_Eg_co_z[i], 
            TRANS_v_In_phys_x[i], TRANS_v_In_phys_y[i], TRANS_v_In_phys_z[i],
            TRANS_cam_phys_x[i], TRANS_cam_phys_y[i], TRANS_cam_phys_z[i],
            TRANS_cam_phys_offset_x[i], TRANS_cam_phys_offset_y[i], TRANS_cam_phys_offset_z[i],
            TRANS_fov_kpc[i], TRANS_box_redshift[i], TRANS_cylrad[i]  ) ; mogfile.write(line) #; print line

    
    mogfile.close()

    
    return 1












if __name__=="__main__":

    moglist = range(305)

    for mog in moglist:
        mogfile = 'hudfmovie_'+str(mog)+'.txt'
        print mogfile
        result = transmogrify('hudfmovie_0.txt',mog,mogfile)


