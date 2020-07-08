#import ez_galaxy
import math
import string
import sys
import struct
import matplotlib
import matplotlib.pyplot as pyplot
#import utils
import numpy as np
#import array
#import astLib.astStats as astStats
import cPickle
import asciitable
import scipy.ndimage
import scipy.stats as ss
import scipy as sp
import pyfits
import make_color_image
import os


# The camera class performs the basic operations.
# It takes as input 10 parameters from the CAMERAX-PARAMETERS HDUs created by Sunrise
# The position and FOV units are in KPC
# It returns an object containing these data plus methods for converting generic
#   x,y,z coordinates from the simulation frame (in Physical kpc!!) into a camera-based coordinate system.
# The camera coordinates are defined such that the axis ranges are [-1,1].
#     The return coordinates can be modified to use a pixel-based grid instead, but this more generic function can be used for both the 
#     CANDELized and perfect images (i.e., on the same axis extent, in matplotlib terms)
# There is one remaining uncertainty -- the sense of the rows and columns in the stored images (and how they are read into a given language).
#     In Python:pyfits/astropy, given simulation coordinates, the returned pixel values correspond to the location on the image given the following assumptions:
#     The "imshow" command was run with origin='lower' and extent=(-1,1,-1,1)
#     The images returned by pyfits from the broadband.fits or _candelized_noise.fits must be **TRANSPOSED** first
#     Presumably there are other iterations of these two settings (or other ways to manipulate the images) that will be satisfactory (or better, if the reason they work is known).
#      -Greg Snyder, 8/21/2014

class camera:
    def __init__(self,x,y,z,dirx,diry,dirz,upx,upy,upz,fov,shape):
        self.x=x
        self.y=y
        self.z=z
        self.dirx=dirx  #x3 unit vector w/ ref to lab frame
        self.diry=diry
        self.dirz=dirz
        self.upx=upx  #x2 unit vector
        self.upy=upy
        self.upz=upz
        self.fov = fov
        #These vectors are defined following the convention at http://en.wikipedia.org/wiki/Pinhole_camera_model
        self.x3vector = np.asarray([self.dirx,self.diry,self.dirz])
        self.x2vector = np.asarray([self.upx,self.upy,self.upz])
        self.x1vector = np.cross(self.x2vector,self.x3vector)
        self.x1vector = self.x1vector/np.linalg.norm(self.x1vector)
        self.shape = shape


        # This is the heart of the matter.  The idea is to express the object's coordinates in the frame of the camera model defined in __init__.
        # Let the object's position expressed in the original frame be A, and the unit vectors i1, i2, i3 be those along the simulation axes.
        # Let the object's position defined in the camera's reference (without shifting the origin yet) be A'.
        # Then, with linear operators, A' = M A, where M is constructed by taking dot products of the camera's unit vectors i' with the original unit vectors i.
        # When the original frame is standard cartesian coords, this equation reduces to the algebra below.
    def express_in_camcoords(self,x,y,z):
        new_x = x*self.x1vector[0] + y*self.x1vector[1] + z*self.x1vector[2]
        new_y = x*self.x2vector[0] + y*self.x2vector[1] + z*self.x2vector[2]
        new_z = x*self.x3vector[0] + y*self.x3vector[1] + z*self.x3vector[2]
        return np.asarray([new_x,new_y,new_z])

    #Wrapper that reconstructs the Sunrise pinhole camera model, expresses the object's position in the camera frame, and computes its position in the image plane.
    def xyz_to_pixelvals(self,x,y,z):
        camdist = (self.x**2 + self.y**2 + self.z**2)**0.5
        camvec = self.express_in_camcoords(x,y,z)
        #define focal length such that image values span from -1 to 1.
        f = camdist/(0.5*self.fov)
        #See guidance at http://en.wikipedia.org/wiki/Pinhole_camera_model
        y1 = (f/camdist)*camvec[0]*1.0
        y2 = (f/camdist)*camvec[1]*1.0
        return y1,y2
    
    #Wrapper that converts a list of Mandelker clump data into list of pixel coordinates, as defined above.
    def process_clump_data(self, clumpdata):
        x = clumpdata['col4']
        y = clumpdata['col5']
        z = clumpdata['col6']
        N = x.shape[0]
        y1_array = np.ndarray(shape=(N))*0.0
        y2_array = np.ndarray(shape=(N))*0.0
        for i in range(N):
            y1,y2 = self.xyz_to_pixelvals(x[i],y[i],z[i])
            y1_array[i] = y1
            y2_array[i] = y2
            
        return y1_array, y2_array



#set_camobj_from_hdu:  Helper function taking a pyfits HDU, the hdu['CAMERAX-PARAMETERS'], and returning an initialized camera object.    
def set_camobj_from_hdu(camparhdu_1):
    camposx = camparhdu_1.header.get('CAMPOSX')
    camposy = camparhdu_1.header.get('CAMPOSY')
    camposz = camparhdu_1.header.get('CAMPOSZ')
    camdirx = camparhdu_1.header.get('CAMDIRX')
    camdiry = camparhdu_1.header.get('CAMDIRY')
    camdirz = camparhdu_1.header.get('CAMDIRZ')
    camupx = camparhdu_1.header.get('CAMUPX')
    camupy = camparhdu_1.header.get('CAMUPY')
    camupz = camparhdu_1.header.get('CAMUPZ')
    fov_kpc = camparhdu_1.header.get('linear_fov')
    shape = camparhdu_1.data.shape
    camobj_1 = camera(camposx,camposy,camposz,camdirx,camdiry,camdirz,camupx,camupy,camupz,fov_kpc,shape)
    return camobj_1



#The main example script demonstrating the use of the camera object.

if __name__=="__main__":
    #Note: You'll need to edit these paths to get the example working.
    clumpcat = '/Users/gfs/Documents/Professional/HydroART_Morphology/VELA26_clumps/clump_catalogue/Nir_clump_cat_a0.420.txt'
    clumpdat = asciitable.read(clumpcat,data_start=6)  #this call skips the "Bulge" clump because there is no whitespace between the first two entries of that line.
    imagefile = '/Users/gfs/Documents/Professional/HydroART_Morphology/VELA26/VELA26_a0.420_0003014__skipir/broadbandz.fits'

    esclump = clumpdat[0]  #for the example above, this one is the "ex-situ" clump
    print esclump
    galx = esclump['col4'] ; galy = esclump['col5'] ; galz = esclump['col6']
    boxx = esclump['col7'] ; boxy = esclump['col8'] ; boxz = esclump['col9']
    
    print galx, galy, galz
    print boxx, boxy, boxz

    next10 = clumpdat[1:]
    #simply the rest of the example catalog entries


    #Open the broadbandz.fits file corresponding to this snapshot
    bb = pyfits.open(imagefile)
    bb.info()

    #Load in the images of interest, but especially the CAMERAX-PARAMETERS HDU
    camparhdu_1 = bb['CAMERA1-PARAMETERS'] #CAMERA1
    camhdu_1 = bb['CAMERA1-BROADBAND']
    #12,15,17 == F606, F850, F125... just for example
    F606_1 = camhdu_1.data[12,:,:]
    F850_1 = camhdu_1.data[15,:,:]
    F125_1 = camhdu_1.data[17,:,:]
    
    #Create camera object from campar.
    camobj_1 = set_camobj_from_hdu(camparhdu_1)
    #This shows how the object works on a single entry of x,y,z (can be anything).
    y1_1,y2_1 = camobj_1.xyz_to_pixelvals(galx,galy,galz)
    #Return values span -1 to 1.
    #This now does the same as the previous line, but for an entire asciitable full of Mandelker clump catalog lines.
    y1_array_1, y2_array_1 = camobj_1.process_clump_data(next10)

    #repeat for cameras of interest
    camparhdu_0 = bb['CAMERA0-PARAMETERS'] #CAMERA1
    camhdu_0 = bb['CAMERA0-BROADBAND']
    F606_0 = camhdu_0.data[12,:,:]
    F850_0 = camhdu_0.data[15,:,:]
    F125_0 = camhdu_0.data[17,:,:]
    camobj_0 = set_camobj_from_hdu(camparhdu_0)
    y1_0,y2_0 = camobj_0.xyz_to_pixelvals(galx,galy,galz)  #these span -1 to 1
    y1_array_0, y2_array_0 = camobj_0.process_clump_data(next10)
    
    camparhdu_6 = bb['CAMERA6-PARAMETERS'] #CAMERA1
    camhdu_6 = bb['CAMERA6-BROADBAND']
    F606_6 = camhdu_6.data[12,:,:]
    F850_6 = camhdu_6.data[15,:,:]
    F125_6 = camhdu_6.data[17,:,:]
    camobj_6 = set_camobj_from_hdu(camparhdu_6)
    y1_6,y2_6 = camobj_6.xyz_to_pixelvals(galx,galy,galz)  #these span -1 to 1
    y1_array_6, y2_array_6 = camobj_6.process_clump_data(next10)


    #Grab corresponding CANDELized images
    candels_dir = '/Users/gfs/Documents/Professional/HydroART_Morphology/VELA26/VELA26_a0.420_0003014__skipir/images/'
    c_0_f606 = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA0-BROADBAND_F606W_candelized_noise.fits'))[0].data
    c_0_f850 = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA0-BROADBAND_F850LP_candelized_noise.fits'))[0].data
    c_0_f125 = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA0-BROADBAND_F125W_candelized_noise.fits'))[0].data
    c_1_f606 = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA1-BROADBAND_F606W_candelized_noise.fits'))[0].data
    c_1_f850 = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA1-BROADBAND_F850LP_candelized_noise.fits'))[0].data
    c_1_f125 = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA1-BROADBAND_F125W_candelized_noise.fits'))[0].data
    c_6_f606 = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA6-BROADBAND_F606W_candelized_noise.fits'))[0].data
    c_6_f850 = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA6-BROADBAND_F850LP_candelized_noise.fits'))[0].data
    c_6_f125 = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA6-BROADBAND_F125W_candelized_noise.fits'))[0].data
    #must convert to reasonable flux units (!!)
    b_ZP = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA0-BROADBAND_F606W_candelized_noise.fits'))[0].header.get('AB_ZP')
    g_ZP = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA0-BROADBAND_F850LP_candelized_noise.fits'))[0].header.get('AB_ZP')
    r_ZP = pyfits.open(os.path.join(candels_dir,'VELA26_a0.420_0003014__skipir_CAMERA0-BROADBAND_F125W_candelized_noise.fits'))[0].header.get('AB_ZP')
    b_ZPfact = 10.0**(-0.4*np.float(b_ZP)) *1.0e14
    g_ZPfact = 10.0**(-0.4*np.float(g_ZP)) *1.0e14
    r_ZPfact = 10.0**(-0.4*np.float(r_ZP)) *1.0e14



    #Create example plot.
    fDISK = pyplot.figure(figsize=(12,8))
    pyplot.subplots_adjust(left=0, right=0.999, bottom=0, top=0.999,wspace=0.0,hspace=0.0)
    DISKnum = fDISK.number

    #########
    #THIS PART I JUST PLAYED AROUND WITH UNTIL IT WORKED
    #IMPORTANT
    bSIM_1 = np.transpose(F606_1) ; gSIM_1 = np.transpose(F850_1) ; rSIM_1 = np.transpose(F125_1)
    bSIM_0 = np.transpose(F606_0) ; gSIM_0 = np.transpose(F850_0) ; rSIM_0 = np.transpose(F125_0)
    bSIM_6 = np.transpose(F606_6) ; gSIM_6 = np.transpose(F850_6) ; rSIM_6 = np.transpose(F125_6)

    #multiply by a factor to convert into nu-based flux units from counts
    bHST_1 = np.transpose(c_1_f606)*b_ZPfact ; gHST_1 = np.transpose(c_1_f850)*g_ZPfact ; rHST_1 = np.transpose(c_1_f125)*r_ZPfact
    bHST_0 = np.transpose(c_0_f606)*b_ZPfact ; gHST_0 = np.transpose(c_0_f850)*g_ZPfact ; rHST_0 = np.transpose(c_0_f125)*r_ZPfact
    bHST_6 = np.transpose(c_6_f606)*b_ZPfact ; gHST_6 = np.transpose(c_6_f850)*g_ZPfact ; rHST_6 = np.transpose(c_6_f125)*r_ZPfact
    #IMPORTANT
    #HERE THERE BE TRANPOSES
    #########


    #to get the colors to look reasonable (started from lambda ratios and modify until happy)
    bscale = 1.0 ; gscale = 6.0/10.0 ; rscale = 6.0/12.5

    #Create RGB images in Lupton+ 04 Scheme
    rgbdata_edge_hires = make_color_image.make_interactive(bSIM_1/bscale,gSIM_1/gscale,rSIM_1/rscale,10.0,9.0)
    rgbdata_face_hires = make_color_image.make_interactive(bSIM_0/bscale,gSIM_0/gscale,rSIM_0/rscale,10.0,9.0)
    rgbdata_6_hires = make_color_image.make_interactive(bSIM_6/bscale,gSIM_6/gscale,rSIM_6/rscale,10.0,9.0)

    rgbdata_edge_lores = make_color_image.make_interactive(bHST_1/bscale,gHST_1*gscale,rHST_1*rscale,0.006,4.0)
    rgbdata_face_lores = make_color_image.make_interactive(bHST_0/bscale,gHST_0*gscale,rHST_0*rscale,0.006,4.0)
    rgbdata_6_lores = make_color_image.make_interactive(bHST_6/bscale,gHST_6*gscale,rHST_6*rscale,0.006,4.0)

            
    #plot 'em.
    axiHST = pyplot.axes([0.0,0.0,0.333,0.5],frameon=True,axisbg='black')
    axiHST.set_xticks([]) ; axiHST.set_yticks([])
    axiHST.imshow(rgbdata_edge_hires,interpolation='nearest',origin='lower',extent=[-1,1,-1,1])  
    axiHST.plot([y1_1],[y2_1],'ok',markerfacecolor='None',markeredgecolor='Lime',markersize=10,alpha=0.75)
    axiHST.plot(y1_array_1,y2_array_1,'ok',markerfacecolor='None',markeredgecolor='Yellow',markersize=5,alpha=0.75)

    axiHST = pyplot.axes([0.333,0.0,0.333,0.5],frameon=True,axisbg='black')
    axiHST.set_xticks([]) ; axiHST.set_yticks([])
    axiHST.imshow(rgbdata_face_hires,interpolation='nearest',origin='lower',extent=[-1,1,-1,1])  
    axiHST.plot([y1_0],[y2_0],'ok',markerfacecolor='None',markeredgecolor='Lime',markersize=10,alpha=0.75)
    axiHST.plot(y1_array_0,y2_array_0,'ok',markerfacecolor='None',markeredgecolor='Yellow',markersize=5,alpha=0.75)

    axiHST = pyplot.axes([0.666,0.0,0.333,0.5],frameon=True,axisbg='black')
    axiHST.set_xticks([]) ; axiHST.set_yticks([])
    axiHST.imshow(rgbdata_6_hires,interpolation='nearest',origin='lower',extent=[-1,1,-1,1])  
    axiHST.plot([y1_6],[y2_6],'ok',markerfacecolor='None',markeredgecolor='Lime',markersize=10,alpha=0.75)
    axiHST.plot(y1_array_6,y2_array_6,'ok',markerfacecolor='None',markeredgecolor='Yellow',markersize=5,alpha=0.75)



    axiHST = pyplot.axes([0.0,0.5,0.333,0.5],frameon=True,axisbg='black')
    axiHST.set_xticks([]) ; axiHST.set_yticks([])
    axiHST.imshow(rgbdata_edge_lores,interpolation='nearest',origin='lower',extent=[-1,1,-1,1])  
    axiHST.plot([y1_1],[y2_1],'ok',markerfacecolor='None',markeredgecolor='Lime',markersize=10,alpha=0.75)
    axiHST.plot(y1_array_1,y2_array_1,'ok',markerfacecolor='None',markeredgecolor='Yellow',markersize=5,alpha=0.75)

    axiHST = pyplot.axes([0.333,0.5,0.333,0.5],frameon=True,axisbg='black')
    axiHST.set_xticks([]) ; axiHST.set_yticks([])
    axiHST.imshow(rgbdata_face_lores,interpolation='nearest',origin='lower',extent=[-1,1,-1,1])  
    axiHST.plot([y1_0],[y2_0],'ok',markerfacecolor='None',markeredgecolor='Lime',markersize=10,alpha=0.75)
    axiHST.plot(y1_array_0,y2_array_0,'ok',markerfacecolor='None',markeredgecolor='Yellow',markersize=5,alpha=0.75)

    axiHST = pyplot.axes([0.666,0.5,0.333,0.5],frameon=True,axisbg='black')
    axiHST.set_xticks([]) ; axiHST.set_yticks([])
    axiHST.imshow(rgbdata_6_lores,interpolation='nearest',origin='lower',extent=[-1,1,-1,1])  
    axiHST.plot([y1_6],[y2_6],'ok',markerfacecolor='None',markeredgecolor='Lime',markersize=10,alpha=0.75)
    axiHST.plot(y1_array_6,y2_array_6,'ok',markerfacecolor='None',markeredgecolor='Yellow',markersize=5,alpha=0.75)

    #save.
    fDISK.savefig('ClumpExample.pdf',format='pdf',dpi=500)
    pyplot.close(fDISK)
