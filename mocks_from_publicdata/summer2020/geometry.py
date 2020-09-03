import math
import string
import sys
import struct
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import array
import scipy.stats as ss
import scipy as sp
#import astropy.io.fits as pyfits
#import cosmocalc as cc
import datetime
#import asciitable
import astropy.io.ascii as ascii
import astropy
import astropy.cosmology


####
#### Name:    geometry.py
#### Author:  Greg Snyder   gsnyder@stsci.edu
#### Purpose:  Generates parameter setup file for mock survey fields from
####           continuous-volume hydrodynamical simulations using Sunrise  (Jonsson 2006).
#### Disclaimer:  This code is provided AS-IS with absolutely NO warranty.
####           It is largely meant as a guide rather than ideal code.
####           It can and should be replaced with other lightcone generation techniques.
####           I make no claims about the immediate usability of this code.
####           That said, I am happy to field questions and discuss issues
####           related to this code. And also to help use it.
#### License:  ?
#### Credit:   Users should cite Snyder et al. (2017)
####           AND ALSO   Kitzbichler & White 2007, Henriques et al. 2013, Overzier et al. 2013, etc.
####           These papers were used to create the lightcone generation algorithm below.
####



class Cosmology:
    def __init__(self, H0=70.0, WM=0.27,WV=0.73,WB=0.0456):
        self.H0=H0
        self.WM=WM
        self.WV=WV
        self.WB=WB
        self.thisc = astropy.cosmology.FlatLambdaCDM(H0=self.H0,Om0=self.WM,Ob0=self.WB)

        self.redshift_grid = np.logspace(-3,2,100)
        self.comoving_mpc_grid=self.thisc.comoving_distance(self.redshift_grid).value
        self.DA_mpc_grid=self.thisc.angular_diameter_distance(self.redshift_grid).value
        #self.comoving_mpc_grid = np.asarray([(cc.cosmocalc(zf,H0=self.H0,WM=self.WM,WV=self.WV))['DCMR_Mpc'] for zf in self.redshift_grid])
        #self.DA_mpc_grid = np.asarray([(cc.cosmocalc(zf,H0=self.H0,WM=self.WM,WV=self.WV))['DA_Mpc'] for zf in self.redshift_grid])

class ReplicatedBox:
    def __init__(self, v_lab, v_ingress):
        self.v_origin=v_lab
        self.v_ingress=v_ingress
        #maybe should define some print/convert functions for this


class LightCone:
    def __init__(self, boxSize, cosmology, name="A Lightcone"):
        self.name=name
        self.cosmology = cosmology
        self.L=boxSize
        self.v1=np.ndarray(shape=(3))
        self.v2=np.ndarray(shape=(3))
        self.v3=np.ndarray(shape=(3))
        self.v4=np.ndarray(shape=(3))
        self.boxlist=[]

    def BasicCone(self, n, m, namelist, zlist, manual_dist_limit=0.0, manual_fov_arcmin=0.0):
        self.n = 1.0*n
        self.m = 1.0*m
        self.namelist = namelist
        self.zlist = zlist
        self.dist_firstRep = np.linalg.norm(np.asarray([self.n,self.m,self.n*self.m])*self.L)
        self.dist_limit = manual_dist_limit
        if manual_dist_limit==0.0:
            self.dist_limit = self.dist_firstRep

        self.redshift_firstRep = np.interp(self.dist_firstRep,self.cosmology.comoving_mpc_grid,self.cosmology.redshift_grid)
        self.numRep = int(self.n*self.m)
        self.x_com = np.asarray( [(self.n - 0.5/self.m)*self.L, (self.n + 0.5/self.m)*self.L] )
        self.y_com = np.asarray( [(self.m - 0.5/self.n)*self.L, (self.m + 0.5/self.n)*self.L] )
        self.z_com = np.asarray([self.n*self.m*self.L])

        self.delta_a_rad = (1.0/(self.n*self.m**2))  #small angle approx?
        self.delta_b_rad = (1.0/(self.m*self.n**2))
        print("WARNING: I'm pretty sure you are assuming that the survey area is small, because I am making some small-angle approximations!  If you are looking for surveys of bigger than ~degree scales, please fix me!")

        self.square_fov_rad = (manual_fov_arcmin/60.0)*(math.pi/180.0)
        if manual_fov_arcmin==0.0:
            self.square_fov_rad = self.delta_b_rad

        self.v1 = np.asarray([(self.x_com)[0],(self.y_com[0]),(self.z_com)[0]])
        self.v2 = np.asarray([(self.x_com)[1],(self.y_com[0]),(self.z_com)[0]])
        self.v3 = np.asarray([(self.x_com)[1],(self.y_com[1]),(self.z_com)[0]])
        self.v4 = np.asarray([(self.x_com)[0],(self.y_com[1]),(self.z_com)[0]])


        self.xaxis = np.asarray([1.0,0.0,0.0])
        self.u3 = np.asarray([(self.n),(self.m),(self.n*self.m)])#/(self.n**2 + self.m**2 + (self.n*self.m)**2)**(0.5)
        self.u3 = self.u3/np.linalg.norm(self.u3)

        self.primaryaxis = np.asarray([0.0,0.0,1.0])
        self.u1 = np.cross(self.u3,self.xaxis)#np.cross(self.xaxis,self.u3)
        self.u1 = self.u1/np.linalg.norm(self.u1)

        self.u2 = np.cross(self.u3,self.u1)
        self.u2 = self.u2/np.linalg.norm(self.u2)

        self.origin=np.asarray([0.0,0.0,0.0])
        self.snapindex = np.where(self.zlist == np.min(self.zlist))

        self.BasicInfo()
        self.ComputeBoxes()

    def BasicInfo(self):
        print("\n")
        print("Information about: ", self.name)
        print("\t Comoving Single Box L = ", self.L)
        print("\t Basic info: n,m = ", self.n, self.m)
        print("\t Approx. Comoving distance at first repeat: ", round(self.dist_firstRep,2))
        print("\t Approx. Redshift at first repeat: ", round(self.redshift_firstRep,2))
        print("\t Number of replications: ", self.numRep)
        print(" ")
        print("\t X range [Mpc] = ", self.x_com)
        print("\t Y range [Mpc] = ", self.y_com)
        print("\t Z height [Mpc] = ", self.z_com)

        print("\n\t del A, arcmin: {:5.2f}".format(self.delta_a_rad*(180.0/math.pi)*60.0))
        print("\t del B, arcmin: {:5.2f}".format(self.delta_b_rad*(180.0/math.pi)*60.0))

        print("\n\t Direction Unit Vector: ", self.u3)
        print("\t Alpha Unit Vector: ", self.u1)
        print("\t Delta Unit Vector: ", self.u2)
        print("\t Test, should be Direction vector: ", np.cross(self.u1,self.u2))

        print(" ")


    def export_runparams(self, filename,follow=False, follow_index=60, swapxy=False , swapxz=False ):

        dirvector = 1.0*self.u3
        alpha_vector = 1.0*self.u1
        delta_vector = 1.0*self.u2
        xind=0
        yind=1
        zind=2

        if swapxy==True:
            temp=dirvector[0]
            dirvector[0]=dirvector[1] ; dirvector[1]=temp
            temp=alpha_vector[0]
            alpha_vector[0]=alpha_vector[1] ; alpha_vector[1]=temp
            temp=delta_vector[0]
            delta_vector[0]=delta_vector[1] ; delta_vector[1]=temp
            xind= 1 ; yind=0 ; zind=2

        if swapxz==True:
            temp=dirvector[0]
            dirvector[0]=dirvector[2] ; dirvector[2]=temp
            temp=alpha_vector[0]
            alpha_vector[0]=alpha_vector[2] ; alpha_vector[2]=temp
            temp=delta_vector[0]
            delta_vector[0]=delta_vector[2] ; delta_vector[2]=temp
            xind= 2 ; yind=1 ; zind=0




        f = open(filename,'w')
        line = '## ' + self.name + ',   LightCone Created, '+ str(datetime.date.today()) + '\n' ; f.write(line) ; print(line)
        line = "## Comoving Single Box L = " + str(self.L) +'\n' ; f.write(line) ; print(line)
        line = "## HubbleParam = " + str(self.cosmology.H0/100.0) + '\n' ; f.write(line) ; print(line) ; h = self.cosmology.H0/100.0
        line = "## Basic info: n,m = " +str( self.n) + " , " + str( self.m) + '\n' ; f.write(line) ; print(line)
        line = "## Approx. Comoving distance at first repeat: " + str( round(self.dist_firstRep,6) ) + '\n' ; f.write(line) ; print(line)
        line = "## Approx. Redshift at first repeat: "  + str( round(self.redshift_firstRep,6) ) + '\n' ; f.write(line) ; print(line)
        line = "## Number of replications: " + str( self.numRep) + '\n' ; f.write(line) ; print(line)
        line = "## del A, arcmin: {:10.5f}".format(self.delta_a_rad*(180.0/math.pi)*60.0) + '\n' ; f.write(line) ; print(line)
        line = "## del B, arcmin: {:10.5f}".format(self.delta_b_rad*(180.0/math.pi)*60.0) + '\n' ; f.write(line) ; print(line)
        line = "## At 0.04 arcsec/pixel, need > {:6.1f} pixels\n".format(self.square_fov_rad*(180.0/math.pi)*3600.0/0.04) ; f.write(line) ; print(line)
        line = "## Direction Unit Vector: " + str( dirvector ) + '\n' ; f.write(line) ; print(line)
        line = "## Alpha Unit Vector: " + str( alpha_vector ) + '\n'  ; f.write(line) ; print(line)
        line = "## Delta Unit Vector: " + str( delta_vector ) + '\n' ; f.write(line) ; print(line)
        line = "## Buffered Cylindricial Radius Maximum: "+str( ((self.boxlist)[-2]).cylinder_radius_approx) + '\n' ; f.write(line) ; print(line)
        line = "## Column 1:  ID#\n" ; f.write(line)
        line = "## Column 2:  Snapshot Label\n" ; f.write(line)
        line = "## Column 3:  Snapshot Redshift\n" ; f.write(line)

        line = "## Column 4:  v_Ingress along x [Comoving h^-1 kpc]\n" ; f.write(line)
        line = "## Column 5:  v_Ingress along y [Comoving h^-1 kpc]\n" ; f.write(line)
        line = "## Column 6:  v_Ingress along z [Comoving h^-1 kpc]\n" ; f.write(line)

        line = "## Column 7:  v_Egress along x [Comoving h^-1 kpc]\n" ; f.write(line)
        line = "## Column 8:  v_Egress along y [Comoving h^-1 kpc]\n" ; f.write(line)
        line = "## Column 9:  v_Egress along z [Comoving h^-1 kpc]\n" ; f.write(line)

        line = "## Column 10:  v_Ingress along x [Physical kpc]\n" ; f.write(line)
        line = "## Column 11:  v_Ingress along y [Physical kpc]\n" ; f.write(line)
        line = "## Column 12:  v_Ingress along z [Physical kpc]\n" ; f.write(line)

        line = "## Column 13:  v_Camera along x [Physical kpc] \n" ; f.write(line)
        line = "## Column 14:  v_Camera along y [Physical kpc] \n" ; f.write(line)
        line = "## Column 15:  v_Camera along z [Physical kpc] \n" ; f.write(line)

        line = "## Column 16:  v_Camera - v_Ingress along x [Physical kpc] \n" ; f.write(line)
        line = "## Column 17:  v_Camera - v_Ingress along y [Physical kpc] \n" ; f.write(line)
        line = "## Column 18:  v_Camera - v_Ingress along z [Physical kpc] \n" ; f.write(line)

        line = "## Column 19:  Square Field of View (smaller axis) at v_Ingress [Physical kpc]\n" ; f.write(line)

        line = "## Column 20:  Geometrically-appropriate redshift at center of box\n" ; f.write(line)
        line = "## Column 21:  Radius buffered to subtend FOV [Comoving h^-1 kpc]\n" ; f.write(line)

        i=0
        MaxRadSize =  ((self.boxlist)[-2]).cylinder_radius_approx
        for box in (self.boxlist)[:-1]:
            if follow==True:
                followbox = (self.boxlist)[follow_index]
            if follow==False:
                followbox=box

            v_in_snap = followbox.v_ingress_local*1000.0*h#np.mod(box.v_ingress, self.L)*1000.0*h  #in comoving kpc h^-1 units
            v_out_snap = followbox.v_egress_local*1000.0*h #np.mod(box.v_egress, self.L)*1000.0*h
            v_in_phys = followbox.v_ingress_local*1000.0/(1.0 + box.mid_z) # in physical kpc
            v_out_phys = followbox.v_egress_local*1000.0/(1.0 + box.mid_z) # in physical kpc

            v_cam_phys = v_in_phys - 1.0*box.camera_offset*1000.0*self.u3/(1.0 + box.mid_z) # in physical kpc, laboratory frame -- does Sunrise translate camera coords too?!?!
            v_cam_cent_phys = v_cam_phys - v_in_phys # in case we want to center on the ingress point
            fov_phys = 2.0*(box.start_distance)*math.sin(self.square_fov_rad/2.0)*1000.0/(1.0 + box.mid_z)  #in physical kpc

            RadSize_snap = box.cylinder_radius_approx*1000.0*h #MaxRadSize*1000.0*h
            line = "{:5d}   {:4s}   {:7.4f}   {:10.4f}   {:10.4f}   {:10.4f}" \
            "   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}" \
            "   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}   {:10.4f}" \
            "   {:10.4f}   {:7.4f}   {:10.4f}\n".format(i,box.snaplabel,box.snapredshift,
                                                                                          v_in_snap[xind], v_in_snap[yind], v_in_snap[zind],
                                                                                          v_out_snap[xind], v_out_snap[yind], v_out_snap[zind],
                                                                                          v_in_phys[xind], v_in_phys[yind], v_in_phys[zind],
                                                                                          v_cam_phys[xind], v_cam_phys[yind], v_cam_phys[zind],
                                                                                          v_cam_cent_phys[xind], v_cam_cent_phys[yind], v_cam_cent_phys[zind],
                                                                                          fov_phys, box.mid_z, RadSize_snap) ; f.write(line) ; print(line)
            i=i+1
        f.close()

    def ComputeBoxes(self):
        print("\t Computing camera parameters for lightcone: ", self.name)

        distancetraveled=0.0
        ingress_point=self.origin
        ingress_snapindex = self.snapindex
        cmpc_from_z0 = np.interp((self.zlist)[ingress_snapindex], self.cosmology.redshift_grid, self.cosmology.comoving_mpc_grid)

        print("cmpc:  ", cmpc_from_z0)

        self.boxlist.append(ReplicatedBox(self.origin,ingress_point))
        i=0
        Nvec = np.asarray([1.0,1.0,1.0])
        volfrac = 0.0
        while (self.dist_limit - distancetraveled > 1e-10):
            box_i = (self.boxlist)[-1]
            box_i.num = i
            testvec = Nvec*np.asarray([self.L,self.L,self.L])  #boundary to test

            ftest = (testvec - box_i.v_ingress)/self.u3  #propagate to nearest boundary
            factor = np.min(ftest)  #how far til we get one exit?
            ind_exit = np.where((ftest - factor) < 1e-10)  #which axis/es was it?
                                                           #print i, ftest, ind_exit[0]
            box_i.v_ingress_local = box_i.v_ingress - (Nvec - 1.0)*self.L
            box_i.v_egress = box_i.v_ingress + factor*self.u3  #this is where the ray leaves this box
            box_i.v_egress_local = box_i.v_egress - (Nvec-1.0)*self.L


            Nvec[ind_exit[0]] = Nvec[ind_exit[0]] + 1.0  #iterate the boundary along these axes; note generically this could be - 1.0 if using arbitrary start/direction


            olddist = distancetraveled
            distancetraveled = np.linalg.norm(box_i.v_egress)
            mid_dist = olddist + (distancetraveled - olddist)/2.0
            mid_z = np.interp(mid_dist,self.cosmology.comoving_mpc_grid, self.cosmology.redshift_grid)

            box_i.far_z = np.interp(distancetraveled,self.cosmology.comoving_mpc_grid, self.cosmology.redshift_grid)
            box_i.near_z = np.interp(olddist,self.cosmology.comoving_mpc_grid, self.cosmology.redshift_grid)
            box_i.mid_z = mid_z  #this is used later
            box_i.mid_dist = mid_dist
            diffs = np.abs(self.zlist - mid_z)
            closest_ind = np.where(diffs == np.min(diffs))   # is this {snapshot selection} the only thing z is used for here?
            box_i.snaplabel = ((self.namelist)[closest_ind[0]])[0]
            box_i.snapredshift = ((self.zlist)[closest_ind[0]])[0]
            box_i.tot_distance_traveled_through = distancetraveled
            box_i.box_distance  = (distancetraveled - olddist)
            box_i.start_distance = olddist
            box_i.camera_offset = box_i.start_distance#/(1.0 + box_i.mid_z) actually, let's keep this in co-moving units #distancetraveled/((1.0 + box_i.snapredshift)) - box_i.box_distance/(1.0 + box_i.snapredshift)  #=~ olddist/(1+z) ...
            box_i.cylinder_radius_approx = ((self.square_fov_rad/2.0)*(2.0**0.5)*1.01)*distancetraveled

            box_i.tot_fov_comoving = (self.square_fov_rad)*distancetraveled  #small angle approx...
            #print closest_ind[0]
            self.boxlist.append(ReplicatedBox((Nvec-1.0)*self.L,box_i.v_egress))  #add the new box
                                                                                  #can update/save some of its basic properties after this

            '''print i, "{:10.3f},  {:10.3f},  {:10.3f},  {:10.3f},  {:12.8f},  {:10.3f},  {:10.3f},  {:5s}".format( np.round_(distancetraveled,3),
                                                                            np.round_(self.delta_b_rad*distancetraveled,3),
                                                                            np.round_(self.delta_b_rad*np.interp(distancetraveled,self.cosmology.comoving_mpc_grid,self.cosmology.DA_mpc_grid),3),np.round_(np.interp(distancetraveled,self.cosmology.comoving_mpc_grid, self.cosmology.redshift_grid), 3), (self.L*np.round_(self.delta_b_rad*distancetraveled,3)**2)/(self.L**3),mid_dist, mid_z, (self.namelist)[closest_ind[0]])'''
            box_i.approx_volume_comoving = (self.L*np.round_(self.square_fov_rad*distancetraveled,3)**2)/(self.L**3)
            volfrac = volfrac + (self.L*np.round_(self.square_fov_rad*distancetraveled,3)**2)/(self.L**3)
            #, np.round_(box_i.v_ingress-box_i.v_origin,3)
            i=i+1

        self.volfrac = volfrac
        #print "Rough Cumulative Volume Fraction (of single box): ", self.volfrac


if __name__=="__main__":

    print("Exploring some things about setting up lightcones...")
    h=0.6774
    L = 20.0/h
    #print "L = ", L, " Mpc"

    #default HUDF-ish lightcone
    n = 15.0 ; m = 14.0
    #print "n,m = ", n,",", m

    fakez = np.logspace(-3,2,100)
    #    comds = np.asarray([(cc.cosmocalc(zf))['DCMR_Mpc'] for zf in fakez])


    delta_a_rad = (1.0/(n*m**2))
    delta_b_rad = (1.0/(m*n**2))

    skyPixel_arcsec = 0.04 #arcsec
    print("ideal ACS-ish scale: {:8.2f}".format(skyPixel_arcsec))
    Npix_A = (delta_a_rad*(180.0/math.pi)*3600.0)/skyPixel_arcsec
    Npix_B = (delta_b_rad*(180.0/math.pi)*3600.0)/skyPixel_arcsec
    print("Npix_A: {:10.1f}".format(Npix_A))
    print("Npix_B: {:10.1f}".format(Npix_B))

    GB_per_slice = 4.0*Npix_A*Npix_B/1e9
    print("GigaBytes per float: {:7.2f}".format(GB_per_slice))

    redshift = np.logspace(-3, 1, 40)
    #print z

    Nz = (redshift.shape)[0]

    #for zi in redshift:
	#    res = cc.cosmocalc(zi)
	#    print "At z= {:6.3f}, D_com= {:6.1f}; DA= {:6.1f} Mpc; DL= {:8.1f}; PS= {:3.1f} kpc/arcsec; dXz= {:5.2f}; dYz= {:5.2f}".format(
	#    round(zi,3), res['DCMR_Mpc'], res['DA_Mpc'], res['DL_Mpc'], res['PS_kpc'], delta_a_rad*res['DCMR_Mpc'], delta_b_rad*res['DCMR_Mpc'])
	#    test = cc.cosmocalc(2.0,H0=71.0,WM=0.27,WV=None)



    #data = asciitable.read('gfm_snaps.txt')
    #data = asciitable.read('snap_v_redshift.txt')

    data = ascii.read('tng_snaps_v_redshift.txt')
    zlist=data['col2'].data
    namelist=np.asarray(data['col1'].data,dtype=str)

    #zlist = np.array(map(float,(data['col2'])))

    #zlist = np.array(map(float,(data['col2'])[:315]))
    #namelist = ((data['col1'])[:315])
    #namelist = np.asarray([(s)[85:89] for s in namelist])

    #namelist = ((data['col1']))
    #namelist = np.asarray([(str(s)) for s in namelist])

    #namelist = np.asarray([(s)[84:88] for s in namelist])

    #print namelist#, zlist
    cosmology = Cosmology(H0=67.74,WM=0.3089,WV=1.0-0.3089,WB=0.0486)
    #hudf_default = LightCone(75.0/h,cosmology,"Default Deep")
    #hudf_default.BasicCone(11.0, 9.0, namelist, zlist)
    #hudf_shallow = LightCone(75.0/h,cosmology,"Default Shallow")
    #hudf_shallow.BasicCone(5.0, 4.0, namelist, zlist)
    #hudf_narrow = LightCone(25.0/h,cosmology,"Default Deep but Narrow")
    #hudf_narrow.BasicCone(11.0, 10.0, namelist, zlist, manual_fov_arcmin=1.0)




    hudf_bigbox_wide = LightCone(75.0/h,cosmology,"Deep 75 Mpc")
    hudf_bigbox_wide.BasicCone(7.0, 6.0, namelist, zlist, manual_dist_limit=10000.0)  #z~8
    hudf_bigbox_wide.export_runparams('tng100_7_6_xyz.txt')
    hudf_bigbox_wide.export_runparams('tng100_7_6_yxz.txt', swapxy=True)
    hudf_bigbox_wide.export_runparams('tng100_7_6_zyx.txt', swapxz=True)


    #hudf_bigbox_vwide = LightCone(75.0/h,cosmology,"Very Wide 75mpc repeated, 136 snaps")
    #hudf_bigbox_vwide.BasicCone(6.0, 5.0, namelist, zlist, manual_dist_limit=11000.0)  #z~18
    #hudf_bigbox_vwide.export_runparams('hudfwide_75Mpc_6_5_xyz.txt')
    #hudf_bigbox_vwide.export_runparams('hudfwide_75Mpc_6_5_yxz.txt', swapxy=True)
    #hudf_bigbox_vwide.export_runparams('hudfwide_75Mpc_6_5_zyx.txt', swapxz=True)


        #hudf_default.export_runparams('hudf_default_75Mpc_11_9_wrongsnaps.txt')



    mpcgrid = cosmology.comoving_mpc_grid
    zgrid = cosmology.redshift_grid

    sizes=[75.0/0.6774,205.0/0.6774,750.0/0.6774,2000.0/0.6774]

    print("{:6s},{:6.0f},{:6.0f},{:6.0f},{:6.0f}".format('box', sizes[0],sizes[1],sizes[2],sizes[3]))
    m = 10.0 ; n = 11.0
    print("{:6.1f},{:6.2f},{:6.2f},{:6.2f},{:6.2f}".format(1.0/(m*n**2.0)*(180.0/math.pi)*60.0,
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[0]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[1]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[2]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[3]),mpcgrid,zgrid)))



    m = 8.0 ; n = 9.0
    print("{:6.1f},{:6.2f},{:6.2f},{:6.2f},{:6.2f}".format(1.0/(m*n**2.0)*(180.0/math.pi)*60.0,
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[0]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[1]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[2]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[3]),mpcgrid,zgrid)))


    m = 6.0 ; n = 7.0
    print("{:6.1f},{:6.2f},{:6.2f},{:6.2f},{:6.2f}".format(1.0/(m*n**2.0)*(180.0/math.pi)*60.0,
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[0]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[1]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[2]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[3]),mpcgrid,zgrid)))


    m = 5.0 ; n = 6.0
    print("{:6.1f},{:6.2f},{:6.2f},{:6.2f},{:6.2f}".format(1.0/(m*n**2.0)*(180.0/math.pi)*60.0,
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[0]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[1]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[2]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[3]),mpcgrid,zgrid)))


    m = 2.0 ; n = 3.0
    print("{:6.1f},{:6.2f},{:6.2f},{:6.2f},{:6.2f}".format(1.0/(m*n**2.0)*(180.0/math.pi),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[0]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[1]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[2]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[3]),mpcgrid,zgrid)))

    m = 2.0 ; n = 1.0
    print("{:6.1f},{:6.2f},{:6.2f},{:6.2f},{:6.2f}".format(1.0/(m*n**2.0)*(180.0/math.pi),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[0]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[1]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[2]),mpcgrid,zgrid),
                                                 np.interp(np.linalg.norm(np.asarray([n,m,n*m])*sizes[3]),mpcgrid,zgrid)))
