# ======================================================================================================================
#       Libraries and Modules

# Core modules
import numpy as np
import matplotlib.pyplot as plt
import illustris_python as ilpy
import translate_coordinates as tc #renato's code for camera projections
import tng_api_utils as tau

# Used only sparingly (maybe remove dependencies?)
import os
import astropy.io.ascii as ascii
import astropy
import astropy.io.fits as fits
import astropy.units as u
from astropy.cosmology import WMAP7,z_at_value
import copy

# Constants 
ilh = tau.tngh # Little H (H_0/100) is set to 0.704
illcos = tau.tngcos # Our cosmology is taken from astropy. 
                    # It uses astropy.cosmology.FlatLambdaCDM(H0=70.4,Om0=0.2726,Ob0=0.0456)
#======================================================================================================================








class lightcone_catalog:
# This class holds the 

    def __init__(self,lightconefile,base_dir,mass_limit=(10.0**9.5),sfr_limit=0.0,mag_limit=None):
        lc_data = ascii.read(lightconefile)
        print("Initializing Lightcone File: ", lightconefile)
        print(lc_data)
        self.lightconefile = lightconefile

        self.cylinder_number = np.int32(lc_data['col1'].data)
        self.snapshot_string = lc_data['col2'].data
        self.snapshot_redshift = lc_data['col3'].data

        self.v_Ingress_x_cmh = lc_data['col4'].data
        self.v_Ingress_y_cmh = lc_data['col5'].data
        self.v_Ingress_z_cmh = lc_data['col6'].data

        self.v_Egress_x_cmh = lc_data['col7'].data
        self.v_Egress_y_cmh = lc_data['col8'].data
        self.v_Egress_z_cmh = lc_data['col9'].data

        self.v_Ingress_x_kpc = lc_data['col10'].data
        self.v_Ingress_y_kpc = lc_data['col11'].data
        self.v_Ingress_z_kpc = lc_data['col12'].data

        self.v_Camera_x_kpc = lc_data['col13'].data
        self.v_Camera_y_kpc = lc_data['col14'].data
        self.v_Camera_z_kpc = lc_data['col15'].data

        self.v_Offset_x_kpc = lc_data['col16'].data
        self.v_Offset_y_kpc = lc_data['col17'].data
        self.v_Offset_z_kpc = lc_data['col18'].data

        self.fov_kpc = lc_data['col19'].data
        self.center_redshift = lc_data['col20'].data
        self.radius_buffer_cmh = lc_data['col21'].data

        xs = None
        xd = None
        self.L_comoving = None

        lines = open(lightconefile,'r')
        for l in lines:
            if "Comoving Single Box L" in l:
                self.L_comoving = np.float32(l.split()[-1])
                self.L_comovingh = round(self.L_comoving*ilh,4)

            if "Delta Unit Vector" in l:
                ss = l.split("[")[-1].split("]")[0].split()
                xs = ss[0]
                ys = ss[1]
                zs = ss[2]

            if "Direction Unit Vector" in l:
                ss = l.split("[")[-1].split("]")[0].split()
                xd = ss[0]
                yd = ss[1]
                zd = ss[2]
            if "del B" in l:
                self.delb_arcmin = np.float32(l.split()[-1])
            if "del A" in l:
                self.dela_arcmin = np.float32(l.split()[-1])
        lines.close()
        assert xs is not None
        assert xd is not None
        assert self.L_comoving is not None

        #camdir
        #just the direction unit vector from the lightcone file
        self.camdir_x = np.float32(xd)
        self.camdir_y = np.float32(yd)
        self.camdir_z = np.float32(zd)

        #camup
        #just the delta unit vector from the lightcone file
        self.camup_x = np.float32(xs)
        self.camup_y = np.float32(ys)
        self.camup_z = np.float32(zs)

        print("    Direction vector: ", self.camdir_x, self.camdir_y, self.camdir_z)
        print("    Up vector: ", self.camup_x, self.camup_y, self.camup_z)
        print("    B FOV, arcmin: ", self.delb_arcmin)
        print("    A FOV, arcmin: ", self.dela_arcmin)
        print("    Ls, Mpc: ", self.L_comoving, self.L_comovingh)

        self.norm_degrees = self.delb_arcmin/60.0

        self.cylinder_object_list = []

        self.base_dir = base_dir

        self.mass_limit = mass_limit
        self.sfr_limit = sfr_limit
        self.mag_limit = mag_limit

        return

#0 - gas
#1 - DM
#4 - stars & WIND
#5 - BHs

    def process_lightcone(self,minz=0.0,maxz=20.0):

        cmd_total = 0.0
        cmx = 0.0
        cmy = 0.0
        cmz = 0.0



        for i,cyl in enumerate(self.cylinder_number):

            cmd_thiscyl = ( (self.v_Egress_x_cmh[i]/ilh - self.v_Ingress_x_cmh[i]/ilh)**2 + (self.v_Egress_y_cmh[i]/ilh - self.v_Ingress_y_cmh[i]/ilh)**2 + (self.v_Egress_z_cmh[i]/ilh - self.v_Ingress_z_cmh[i]/ilh)**2 )**0.5
            cmd_begin = cmd_total
            cmd_end = cmd_begin + cmd_thiscyl
            cmd_total = cmd_end

            cz=self.center_redshift[i]


            #world coordinates of ingress point
            cmx_begin = 1.0*cmx
            cmy_begin = 1.0*cmy
            cmz_begin = 1.0*cmz

            #world coordinates of egress points
            cmx = cmx_begin + (self.v_Egress_x_cmh[i]/ilh - self.v_Ingress_x_cmh[i]/ilh)
            cmy = cmy_begin + (self.v_Egress_y_cmh[i]/ilh - self.v_Ingress_y_cmh[i]/ilh)
            cmz = cmz_begin + (self.v_Egress_z_cmh[i]/ilh - self.v_Ingress_z_cmh[i]/ilh)

            if i > 1000:
                continue

            if cz < minz:
                continue
            if cz > maxz:
                continue


            testf = 'test_'+str(cyl)+'.pdf'
            f1 = plt.figure(figsize=(10.5,10.5), dpi=150)
            plt.subplots_adjust(left=0.11, right=0.98, bottom=0.08, top=0.99,wspace=0.25,hspace=0.25)
            skip = 500

            #determine snapshot of interest
            print("Processing Cylinder: ", cyl, i, self.snapshot_redshift[i])
            snapnum = self.snapshot_string[i]

            #old corrupt snaps for illustris-1
            #if snapnum==53:
            #    snapnum=52
            #if snapnum==55:
            #    snapnum=54

            print("    Snapshot Number: ", snapnum)

            #load subhalo catalogs for this snapshot

            fields=['SubhaloMass','SubhaloMassInMaxRad','SubhaloMassInRadType','SubhaloMassInMaxRadType','SubhaloPos','SubhaloSFR','SubhaloSFRinRad','SubhaloVel','SubhaloBHMass','SubhaloBHMdot','SubhaloStellarPhotometrics','SubhaloWindMass']
            subhalos = ilpy.groupcat.loadSubhalos(self.base_dir,snapnum,fields=fields)
            print("    Loaded subhalos: ", subhalos['count'], subhalos['SubhaloMassInRadType'].shape)


            subhalos = self.periodicize(subhalos,self.L_comovingh*1000.0)
            
            mstar_msun = subhalos['SubhaloMassInRadType'][:,4]*(1.0e10)/ilh
            mgas_msun = subhalos['SubhaloMassInRadType'][:,0]*(1.0e10)/ilh #includes wind mass
            mbh_msun = subhalos['SubhaloMassInRadType'][:,5]*(1.0e10)/ilh

            baryonmass_msun = mstar_msun + mgas_msun + mbh_msun #within 2x stellar half mass radius... best?

            mhalo_msun = subhalos['SubhaloMass']*(1.0e10)/ilh

            sfr = subhalos['SubhaloSFR']*1.0

            gmag_ABabs=subhalos['SubhaloStellarPhotometrics'][:,4]*1.0
            distmod=illcos.distmod(cz).value
            gmag=gmag_ABabs+distmod




            #cull liberally by mass... how?
            #here, either massive OR star-forming!  probably must refine this.. how big?

            if self.mag_limit is None:
                mi = np.where(np.logical_and(baryonmass_msun > self.mass_limit, sfr >self.sfr_limit))[0]
            else:
                mi = np.where(np.logical_and(gmag < self.mag_limit,baryonmass_msun > 0.0))[0]


            if mi.shape[0]==0:
                cylinder_obj = None
                self.cylinder_object_list.append(cylinder_obj)
                continue
                f1.close()

            print("    Selected number: ", mi.shape)
            print("    Mstar statistics: ", np.min(mstar_msun[mi]), np.max(mstar_msun[mi]), np.median(mstar_msun[mi]))
            print("    Mgas  statistics: ", np.min(mgas_msun[mi]), np.max(mgas_msun[mi]), np.median(mgas_msun[mi]))
            print("    Mag  statistics : ", np.min(gmag[mi]), np.max(gmag[mi]), np.median(gmag[mi]))



            #mi is index into subhalo catalog -- the critical index numbers we need to save!!!
            xpos = subhalos['SubhaloPos'][mi,0] #in cKpc/h of max bound part
            ypos = subhalos['SubhaloPos'][mi,1]
            zpos = subhalos['SubhaloPos'][mi,2]


            #project geometry

            #campos
            #in phys kpc, offset values from lightcone file!
            xoff = self.v_Offset_x_kpc[i]
            yoff = self.v_Offset_y_kpc[i]
            zoff = self.v_Offset_z_kpc[i]

            #the position here I think doesn't matter???
            camera = tc.Camera([0,0,0],[self.camdir_x,self.camdir_y,self.camdir_z],[self.camup_x,self.camup_y,self.camup_z])

            #galaxy world position
            #convert to phys kpc following Renato's lead in translate_coordinates.py
            #note there's an extra translation in the sunrise calcs, so we can discard that here

            #box coordinates relative to ingress coordinate
            boxX = (xpos/ilh) - self.v_Ingress_x_cmh[i]/ilh
            boxY = (ypos/ilh) - self.v_Ingress_y_cmh[i]/ilh
            boxZ = (zpos/ilh) - self.v_Ingress_z_cmh[i]/ilh

            axi = f1.add_subplot(2,2,1)
            axi.set_ylabel('boxI X',size=7,labelpad=1)
            axi.set_xlabel('boxI Z',size=7,labelpad=1)
            axi.tick_params(axis='both',which='major',labelsize=7)
            axi.plot(boxZ[::skip],boxX[::skip],'ok')

            #add box coordinate to world coordinate of ingress point
            worldX = boxX+cmx_begin
            worldY = boxY+cmy_begin
            worldZ = boxZ+cmz_begin

            axi = f1.add_subplot(2,2,2)
            axi.set_ylabel('world X',size=7,labelpad=1)
            axi.set_xlabel('world Z',size=7,labelpad=1)
            axi.tick_params(axis='both',which='major',labelsize=7)
            axi.plot(worldZ[::skip],worldX[::skip],'ok')
            axi.plot([np.min(worldZ),np.max(worldZ)],[cmx_begin,cmx_begin],color='red')



            velX = subhalos['SubhaloVel'][mi,0]
            velY = subhalos['SubhaloVel'][mi,1]
            velZ = subhalos['SubhaloVel'][mi,2]

            #galaxy cam position, in comoving kpc
            galaxy_camera_posx,galaxy_camera_posy,galaxy_camera_posz = camera.cameraCoordinates_vector(worldX,worldY,worldZ)
            galaxy_camera_velx,galaxy_camera_vely,galaxy_camera_velz = camera.cameraCoordinates_vector(velX,velY,velZ)

            axi = f1.add_subplot(2,2,3)
            axi.set_ylabel('cam X',size=7,labelpad=1)
            axi.set_xlabel('cam Z',size=7,labelpad=1)
            axi.tick_params(axis='both',which='major',labelsize=7)
            axi.plot(galaxy_camera_posz[::skip],galaxy_camera_posx[::skip],'ok')


            #galaxy projection using spherical coords
            y1 = np.arctan2(galaxy_camera_posx,galaxy_camera_posz)/(0.5*(self.delb_arcmin/60.0)*(np.pi/180.0))
            y2 = np.arctan2(galaxy_camera_posy,galaxy_camera_posz)/(0.5*(self.delb_arcmin/60.0)*(np.pi/180.0))
            #range = [-1,1] = FOV = self.norm_degrees



            axi = f1.add_subplot(2,2,4)
            axi.set_ylabel('cam Y1',size=7,labelpad=1)
            axi.set_xlabel('cam Y2',size=7,labelpad=1)
            axi.set_xlim(-3,3)
            axi.set_ylim(-3,3)
            axi.tick_params(axis='both',which='major',labelsize=7)
            axi.plot(y1,y2,'ok',markersize=0.5,mew=0.0)
            axi.plot([-1,-1],[-1,1],color='red')
            axi.plot([-1,1],[-1,-1],color='red')
            axi.plot([1,1],[-1,1],color='red')
            axi.plot([-1,1],[1,1],color='red')


            #all values correspond to mi vector
            #cull by RA, DEC, and segment length
            ci = np.where(np.logical_and(np.logical_and(np.logical_and(np.abs(y1) <= 1.0, np.abs(y2) <= 1.0),galaxy_camera_posz <= cmd_end),galaxy_camera_posz > cmd_begin))[0]
            print("    Selected N galaxies in FOV: ", ci.shape)
            axi.plot(y1[ci],y2[ci],'or',markersize=0.7,mew=0.0)



            RA_deg = y1[ci]*self.norm_degrees/2.0
            DEC_deg = y2[ci]*self.norm_degrees/2.0

            #save interesting quantities
            if ci.shape[0] > 0:
                print(cyl, cmd_begin, np.min(galaxy_camera_posz[ci]))
                print(cyl, cmd_end, np.max(galaxy_camera_posz[ci]))

                cylinder_obj = cylinder_catalog(snapnum,subhalos,mi,ci,RA_deg,DEC_deg, self.snapshot_redshift[i],
                                                galaxy_camera_posx,galaxy_camera_posy,galaxy_camera_posz,self.center_redshift[i],
                                                galaxy_camera_velx,galaxy_camera_vely,galaxy_camera_velz,cyl,gmag[mi])
            else:
                cylinder_obj = None

            self.cylinder_object_list.append(cylinder_obj)

            #f1.savefig(testf)
            plt.close(f1)

        return self

    def periodicize(self,subhalos,boxL):
        xpos = subhalos['SubhaloPos'][:,0].flatten() #in cKpc/h of max bound part
        ypos = subhalos['SubhaloPos'][:,1].flatten()
        zpos = subhalos['SubhaloPos'][:,2].flatten()

        N = xpos.shape[0]

        sid = np.arange(N)
        
        new_subhalos = copy.copy(subhalos)

        new_x = copy.copy(xpos)
        new_y = copy.copy(ypos)
        new_z = copy.copy(zpos)


        new_subhalos['SubFindID'] = np.concatenate((sid,sid))
        new_subhalos['SubFindID'] = np.concatenate((new_subhalos['SubFindID'],sid))
        new_subhalos['SubFindID'] = np.concatenate((new_subhalos['SubFindID'],sid))
        new_subhalos['SubFindID'] = np.concatenate((new_subhalos['SubFindID'],sid))
        new_subhalos['SubFindID'] = np.concatenate((new_subhalos['SubFindID'],sid))
        new_subhalos['SubFindID'] = np.concatenate((new_subhalos['SubFindID'],sid))

        
        keys = subhalos.keys()
        for key in keys:
            if key=='SubhaloPos':
                #special
                #x repeat

                new_x = np.concatenate((new_x,xpos+boxL))
                new_x = np.concatenate((new_x,xpos-boxL))
                new_y = np.concatenate((new_y,ypos))
                new_y = np.concatenate((new_y,ypos))
                new_z = np.concatenate((new_z,zpos))
                new_z = np.concatenate((new_z,zpos))
                #y repeat
                new_x = np.concatenate((new_x,xpos))
                new_x = np.concatenate((new_x,xpos))
                new_y = np.concatenate((new_y,ypos+boxL))
                new_y = np.concatenate((new_y,ypos-boxL))
                new_z = np.concatenate((new_z,zpos))
                new_z = np.concatenate((new_z,zpos))
                #z repeat
                new_x = np.concatenate((new_x,xpos))
                new_x = np.concatenate((new_x,xpos))
                new_y = np.concatenate((new_y,ypos))
                new_y = np.concatenate((new_y,ypos))
                new_z = np.concatenate((new_z,zpos+boxL))
                new_z = np.concatenate((new_z,zpos-boxL))

                new_pos = np.column_stack((new_x,new_y,new_z))

                new_subhalos[key] = new_pos

            elif key=='count':
                new_subhalos[key] = 7*subhalos[key]
            else:
                new_subhalos[key] = np.concatenate((new_subhalos[key],subhalos[key]))
                new_subhalos[key] = np.concatenate((new_subhalos[key],subhalos[key]))
                new_subhalos[key] = np.concatenate((new_subhalos[key],subhalos[key]))
                new_subhalos[key] = np.concatenate((new_subhalos[key],subhalos[key]))
                new_subhalos[key] = np.concatenate((new_subhalos[key],subhalos[key]))
                new_subhalos[key] = np.concatenate((new_subhalos[key],subhalos[key])) #7 total boxes

        return new_subhalos

    def output_catalog(self,outfile):
        print("    Saving catalog: ", outfile)

        fobj = open(outfile,'w')

        fobj.write('## Lightcone Catalog File for input geometry: '+self.lightconefile+'\n')
        fobj.write('## Catalog source directory: '+self.base_dir+'\n')
        fobj.write('## Square FOV (arcmin): {:12.6f}'.format(self.delb_arcmin)+'\n')
        fobj.write('## Area (arcmin^2): {:12.6f}'.format(self.delb_arcmin**2)+'\n')
        fobj.write('## Baryonic Mass Lower Limit (Msun) : {:10.5e}'.format(self.mass_limit)+'\n')
        fobj.write('## Assumed Cosmology: '+WMAP7.__str__()+'\n')
        fobj.write('## Creator:  Teddy Pena (STScI) \n')
        fobj.write('## Catalog & Data Release Reference:  Nelson et al. (2019) \n')
        fobj.write('## Catalog & Data Release URL: tng-project.org/data \n')
        fobj.write('## Column 01: Snapshot number \n')
        fobj.write('## Column 02: Subhalo Index \n')
        fobj.write('## Column 03: RA (degrees) \n')
        fobj.write('## Column 04: DEC (degrees) \n')
        fobj.write('## Column 05: RA (proper kpc at true z) \n')
        fobj.write('## Column 06: DEC (proper kpc at true z) \n')
        fobj.write('## Column 07: RA (proper kpc at inferred z) \n')
        fobj.write('## Column 08: DEC (proper kpc at inferred z) \n')
        fobj.write('## Column 09: True cosmological redshift \n')
        fobj.write('## Column 10: Inferred redshift (includes peculiar v) \n')
        fobj.write('## Column 11: Peculiar redshift; Peculiar Velocity / Speed of Light \n')
        fobj.write('## Column 12: True scale at cosmological z, in kpc/arcsec \n')
        fobj.write('## Column 13: [Mpc] Comoving X in Observer Coordinates \n')
        fobj.write('## Column 14: [Mpc] Comoving Y in Observer Coordinates \n')
        fobj.write('## Column 15: [Mpc] Comoving Z in Observer Coordinates \n')
        fobj.write('## Column 16: [Mpc] True Angular Diameter Distance to observer \n')
        fobj.write('## Column 17: [Mpc] Inferred Angular Diameter Distance to observer \n')
        fobj.write('## Column 18: Snapshot redshift \n')
        fobj.write('## Column 19: Geometrically appropriate redshift at center of this cylinder \n')
        fobj.write('## Column 20: Lightcone cylinder number \n')
        fobj.write('## Column 21: [Msun] Stellar mass within 2X stellar half mass radius\n')
        fobj.write('## Column 22: [Msun] Total gas mass within 2X stellar half mass radius\n')
        fobj.write('## Column 23: [Msun] Total mass of this subhalo (excludes children subhalos) \n')
        fobj.write('## Column 24: [Msun] Total BH mass within 2X stellar half mass radius\n')
        fobj.write('## Column 25: [Msun] Total baryon mass within 2X stellar half mass radius\n')
        fobj.write('## Column 26: [Msun/year] SFR within 2X stellar half mass radius\n')
        fobj.write('## Column 27: [(10^10 Msun/h) / (0.978 Gyr/h)] Total BH accretion rate within subhalo\n')
        fobj.write('## Column 28: [Mpc] Camera X in Observer Coordinates (Proper X at z; a transverse coordinate) \n')
        fobj.write('## Column 29: [Mpc] Camera Y in Observer Coordinates (Proper Y at z; a transverse coordinate)\n')
        fobj.write('## Column 30: [Mpc] Camera Z in Observer Coordinates (Proper Z at z; should be almost exactly Column 16)\n')
        fobj.write('## Column 31: [AB Mag] Intrinsic stellar g absolute magnitude (BC03) \n')
        fobj.write('## Column 32: [AB Mag] Intrinsic stellar r absolute magnitude (BC03) \n')
        fobj.write('## Column 33: [AB Mag] Intrinsic stellar i absolute magnitude (BC03) \n')
        fobj.write('## Column 34: [AB Mag] Intrinsic stellar z absolute magnitude (BC03) \n')
        fobj.write('## Column 35: [km/s] Galaxy motion in transverse Camera X direction \n')
        fobj.write('## Column 36: [km/s] Galaxy motion in transverse Camera Y direction \n')
        fobj.write('## Column 37: [km/s] Galaxy motion in line-of-sight Camera Z direction ; the Peculiar Velocity \n')
        fobj.write('## Column 38: [km/s] Cosmological expansion velocity at true z (Column 10 measures Column 37+38)\n')
        fobj.write('## Column 39: [AB Mag] Apparent total rest-frame g-band magnitude (BC03) \n')

        for cylobj in self.cylinder_object_list:
            if cylobj is not None:
                cylobj.print_cylinder(fobj)

        fobj.close()
        return


class cylinder_catalog:
    def __init__(self,snapnum,subhalos,mi,ci,RA_deg,DEC_deg,snapz,galaxy_camera_posx,galaxy_camera_posy,galaxy_camera_posz,centerz,galaxy_camera_velx,galaxy_camera_vely,galaxy_camera_velz,cyl,gmag):

        #fields=['SubhaloMass','SubhaloMassInMaxRad','SubhaloMassInRadType','SubhaloMassInMaxRadType','SubhaloPos','SubhaloSFR','SubhaloSFRinRad','SubhaloVel','SubhaloBHMass','SubhaloBHMdot','SubhaloStellarPhotometrics','SubhaloWindMass']

        self.snapshot_number = snapnum + np.zeros_like(ci)

        self.subhalo_index = subhalos['SubFindID'][mi[ci]]

        self.RA_deg = RA_deg
        self.DEC_deg = DEC_deg

        self.snapz = snapz + np.zeros_like(RA_deg)
        self.center_z = centerz + np.zeros_like(RA_deg)

        self.cylinder_number = cyl + np.zeros_like(self.subhalo_index)

        self.galaxy_comoving_x_mpc = galaxy_camera_posx[ci]/1000.0
        self.galaxy_comoving_y_mpc = galaxy_camera_posy[ci]/1000.0
        self.galaxy_comoving_z_mpc = galaxy_camera_posz[ci]/1000.0

        #self.galaxy_camera_posx = galaxy_camera_posx[ci]/1000.0
        #self.galaxy_camera_posy = galaxy_camera_posy[ci]/1000.0
        #self.galaxy_camera_posz = galaxy_camera_posz[ci]/1000.0

        self.galaxy_camera_velx = galaxy_camera_velx[ci]
        self.galaxy_camera_vely = galaxy_camera_vely[ci]
        self.galaxy_camera_velz = galaxy_camera_velz[ci]

        self.galaxy_peculiar_vr = 1.0*self.galaxy_camera_velz
        self.galaxy_peculiar_z = 1.0*self.galaxy_peculiar_vr/(astropy.constants.c.value/1.0e3)

        self.cosmological_redshift = np.zeros_like(self.RA_deg)

        for i,index in enumerate(ci):
            self.cosmological_redshift[i] = np.float64(z_at_value(WMAP7.comoving_distance, self.galaxy_comoving_z_mpc[i]*u.megaparsec,ztol=1e-12,maxfun=2000))

        self.hubble_velocity = self.cosmological_redshift*astropy.constants.c.value/1.0e3 #in km/s

        self.galaxy_observed_z = 1.0*self.cosmological_redshift + self.galaxy_peculiar_z

        #self.galaxy_comoving_x_mpc = self.galaxy_camera_posx*(1.0 + self.cosmological_redshift)
        #self.galaxy_comoving_y_mpc = self.galaxy_camera_posy*(1.0 + self.cosmological_redshift)
        #self.galaxy_comoving_z_mpc = self.galaxy_camera_posz*(1.0 + self.cosmological_redshift)

        self.galaxy_camera_posx = self.galaxy_comoving_x_mpc/(1.0 + self.cosmological_redshift)
        self.galaxy_camera_posy = self.galaxy_comoving_y_mpc/(1.0 + self.cosmological_redshift)
        self.galaxy_camera_posz = self.galaxy_comoving_z_mpc/(1.0 + self.cosmological_redshift)

        self.angdiam_mpc = np.asarray(WMAP7.angular_diameter_distance(self.cosmological_redshift))
        self.kpc_per_arcsec = np.asarray(WMAP7.kpc_proper_per_arcmin(self.cosmological_redshift)/60.0)

        self.observed_angdiam_mpc = np.asarray(WMAP7.angular_diameter_distance(self.galaxy_observed_z))
        self.observed_comoving_mpc = np.asarray(WMAP7.comoving_distance(self.galaxy_observed_z))
        self.observed_kpc_per_arcsec = np.asarray(WMAP7.kpc_proper_per_arcmin(self.galaxy_observed_z)/60.0)

        self.RA_kpc = self.RA_deg*3600.0*self.kpc_per_arcsec
        self.DEC_kpc = self.DEC_deg*3600.0*self.kpc_per_arcsec
        self.observed_RA_kpc = self.RA_deg*3600.0*self.observed_kpc_per_arcsec
        self.observed_DEC_kpc = self.DEC_deg*3600.0*self.observed_kpc_per_arcsec

        self.mstar_msun = subhalos['SubhaloMassInRadType'][self.subhalo_index,4]*(1.0e10)/ilh
        self.mgas_msun = subhalos['SubhaloMassInRadType'][self.subhalo_index,0]*(1.0e10)/ilh #includes wind mass
        self.mbh_msun = subhalos['SubhaloMassInRadType'][self.subhalo_index,5]*(1.0e10)/ilh
        self.mhalo_msun = subhalos['SubhaloMass'][self.subhalo_index]*(1.0e10)/ilh

        self.baryonmass_msun = self.mstar_msun + self.mgas_msun + self.mbh_msun #within 2x stellar half mass radius... best?

        self.xpos_ckh = subhalos['SubhaloPos'][self.subhalo_index,0] #in cKpc/h of max bound part
        self.ypos_ckh = subhalos['SubhaloPos'][self.subhalo_index,1]
        self.zpos_ckh = subhalos['SubhaloPos'][self.subhalo_index,2]

        self.xpos_pmpc = (self.xpos_ckh*1.0/(1.0 + snapz )/ilh)/1.0e3
        self.ypos_pmpc = (self.ypos_ckh*1.0/(1.0 + snapz )/ilh)/1.0e3
        self.zpos_pmpc = (self.zpos_ckh*1.0/(1.0 + snapz )/ilh)/1.0e3

        self.xvel_kms = subhalos['SubhaloVel'][self.subhalo_index,0]
        self.yvel_kms = subhalos['SubhaloVel'][self.subhalo_index,1]
        self.zvel_kms = subhalos['SubhaloVel'][self.subhalo_index,2]

        self.sfr = subhalos['SubhaloSFRinRad'][self.subhalo_index]
        self.bhmdot = subhalos['SubhaloBHMdot'][self.subhalo_index]

        self.gmag = subhalos['SubhaloStellarPhotometrics'][self.subhalo_index,4]
        self.rmag = subhalos['SubhaloStellarPhotometrics'][self.subhalo_index,5]
        self.imag = subhalos['SubhaloStellarPhotometrics'][self.subhalo_index,6]
        self.zmag = subhalos['SubhaloStellarPhotometrics'][self.subhalo_index,7]

        self.gmag_apparent=gmag[ci]

        #self.total_redshift =

        return


    def print_cylinder(self,outobj):

        for i,shi in enumerate(self.subhalo_index):

            thisline = '{:8d}{:12d}  {:12.6f}  {:12.6f}  {:10.2f}  {:10.2f}  {:10.2f}  {:10.2f}  '\
                       '{:12.8f}  {:12.8f}  {:12.4e}  {:8.4f}  '\
                       '{:10.4f}  {:10.4f}  {:16.4f}  {:16.4f}  {:16.4f}  {:12.8f}  {:12.8f}  {:8d}'\
                       '{:12.4e}  {:12.4e}  {:12.4e}  {:12.4e}  {:12.4e}  {:16.4f}  {:10.4e}'\
                       '  {:10.4f}  {:10.4f}  {:16.4f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}    {:8.2f}  {:8.2f}  {:8.2f}  {:12.4e}  {:8.2f}'\
                       '\n'.format(self.snapshot_number[i],shi,self.RA_deg[i],self.DEC_deg[i],self.RA_kpc[i],self.DEC_kpc[i],self.observed_RA_kpc[i],self.observed_DEC_kpc[i],
                                   self.cosmological_redshift[i],self.galaxy_observed_z[i],self.galaxy_peculiar_z[i],self.kpc_per_arcsec[i],
                                   self.galaxy_comoving_x_mpc[i],self.galaxy_comoving_y_mpc[i],self.galaxy_comoving_z_mpc[i],self.angdiam_mpc[i],self.observed_angdiam_mpc[i],self.snapz[i],self.center_z[i],self.cylinder_number[i],
                                   self.mstar_msun[i],self.mgas_msun[i],self.mhalo_msun[i],self.mbh_msun[i],self.baryonmass_msun[i],self.sfr[i],self.bhmdot[i],
                                   self.galaxy_camera_posx[i],self.galaxy_camera_posy[i],self.galaxy_camera_posz[i],
                                   self.gmag[i],self.rmag[i],self.imag[i],self.zmag[i],
                                   self.galaxy_camera_velx[i],self.galaxy_camera_vely[i],self.galaxy_camera_velz[i],self.hubble_velocity[i],self.gmag_apparent[i])
            outobj.write(thisline)


        return


def process_lightcone_catalog(lightcone=None,base_dir=None,mass_limit=10.0**9.5,sfr_limit=0.0,mag_limit=None):
    assert (lightcone is not None) and (base_dir is not None)
    assert os.path.lexists(base_dir)

    catalog_object = lightcone_catalog(lightcone,base_dir,mass_limit=mass_limit,sfr_limit=sfr_limit,mag_limit=mag_limit)




    return catalog_object




if __name__=="__main__":


    #start with conservative limits -- should be relatively few sources

    magl=30
    minz=0.1
    maxz=8.8

    #will need to modify inputs below, as well as codes above

    #currently uses local Illustris galaxy catalog data, but I would prefer to use API access calls (if not too slow) or JupyterLab

    #let's start with TNG300-3 simulation for data volume reasons?

    catalog_xyz = process_lightcone_catalog(lightcone="./tng300_6_5_xyz.txt",base_dir='/home/tnguser/sims.TNG/TNG300-1/output/',mag_limit=magl)
    catalog_xyz = catalog_xyz.process_lightcone(minz=minz,maxz=maxz)
    catalog_xyz.output_catalog('./Lightcone_TNG300-1_mag30_6_5_xyz.txt')