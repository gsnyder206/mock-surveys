# ======================================================================================================================
#       Libraries and Modules
from astropy.io import ascii
from astropy.io import fits

import numpy as np
import argparse
import math
# ======================================================================================================================



# ======================================================================================================================
#
#       This script was written at the Space Telescope Science Institute during a Summer Internship Program in June 2015
#       funded by the Institute of International Education (IIE) during the Brazilian Scientific Mobility Program (BSMP)
#
#       Its main function is to perform the translation and rotation of galaxies coordinates in the Illustris Simulation
#       using a Pinhole Camera Model to identify the Galaxies within the simulation and link their properties with the
#       observed galaxies in the final RGB image.
#
#      Author:          Renato da Silva Guimaraes <r.guimaraess@gmail.com>
#      Supervisor:      Gregory F. Snyder <gsnyder@stsci.edu>
#
#
#       Further work, corrections and changes to make the script more accurate were made at the Kavli Institute for
#       Astrophysics and Space Research at the MIT.
#
#       Supervisor:     Paul Torrey <ptorrey@mit.edu>
#
#
#      HOW TO USE
#
#      $ python translate_coordinates.py -config FIELDC_SB27.conf
#
#      For information on how to create the configuration file, please read the manual or the reader in the
#      example configuration file.
#
# ======================================================================================================================


# ======================================================================================================================
#       This class was based on Greg's code, it performs the translation and rotation to the camera reference frame and
#       converts it to pixel coordinates.
# ======================================================================================================================
class Camera():
    def __init__(self, cam_pos, cam_dir, cam_up):

        """
        :param cam_pos: translation vector with the camera position, in cartesian coordinates
        :param cam_dir: unit vector in the direction the camera is poiting, in cartesian coordinates
        :param cam_up: unit vector in the direction considered up for the camera, in cartesian coordinates
        :return: set the variables for inner use and perform the cross product of cam_dir with cam_up to create the third axis
        """

        self.x = cam_pos[0]
        self.y = cam_pos[1]
        self.z = cam_pos[2]
        self.xdir = cam_dir[0]
        self.ydir = cam_dir[1]
        self.zdir = cam_dir[2]
        self.xup = cam_up[0]
        self.yup = cam_up[1]
        self.zup = cam_up[2]

        # Creating the camera base vectors
        self.u3_camera = np.array([self.xdir, self.ydir, self.zdir])
        self.u2_camera = np.array([self.xup, self.yup, self.zup])
        self.u1_camera = np.cross(self.u2_camera, self.u3_camera)
        self.u1_camera = self.u1_camera/np.linalg.norm(self.u1_camera)

    def cameraCoordinates(self, world_vector):
        """
        :param world_vector: Object position in the world reference frame, in cartesian coordinates
        :return: Objects (x, y, z) coordinates translated to the camera position
        """

        worldX = world_vector[0]
        worldY = world_vector[1]
        worldZ = world_vector[2]

        x = worldX*self.u1_camera[0] + worldY*self.u1_camera[1] + worldZ*self.u1_camera[2]
        y = worldX*self.u2_camera[0] + worldY*self.u2_camera[1] + worldZ*self.u2_camera[2]
        z = worldX*self.u3_camera[0] + worldY*self.u3_camera[1] + worldZ*self.u3_camera[2]

        return np.array([x,y,z])

    #Greg vectorized this code on 1/5/2016
    def cameraCoordinates_vector(self, worldX,worldY,worldZ):
        """
        :param world_vector: Object position in the world reference frame, in cartesian coordinates
        :return: Objects (x, y, z) coordinates translated to the camera position
        """

        x = worldX*self.u1_camera[0] + worldY*self.u1_camera[1] + worldZ*self.u1_camera[2]
        y = worldX*self.u2_camera[0] + worldY*self.u2_camera[1] + worldZ*self.u2_camera[2]
        z = worldX*self.u3_camera[0] + worldY*self.u3_camera[1] + worldZ*self.u3_camera[2]

        return x,y,z


    # This is only used if yo want to make the projection from cartesian coordinates, which is not recommended because
    # it introduces more distortion than a spherical projection.
    def pixelCoordinates(self, world_coordinates, fov):
        """
        :param cam_x: Objects X coordinate in the camera reference frame
        :param cam_y: Objects Y coordinate in the camera reference frame
        :param cam_z: Objects Z coordinate in the camera reference frame
        :param fov: Field of view in the current cube
        :return: coordinate in pixels, in a range [-1,1]
        """
        camera_coordinates = self.cameraCoordinates(world_coordinates)

        focal_length = 1.0/(0.5*fov)

        result = np.array([camera_coordinates[0]*focal_length, camera_coordinates[1]*focal_length])
        return result


# ======================================================================================================================
#       Reads the configuration file and sets the default variables if necessary
# ======================================================================================================================
def read_configuration(filename):

    configuration_file = ascii.read(filename)
    error = ''

    if 'OUTPUT1' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('OUTPUT1')
        output1 = configuration_file['VALUE'][index]

    else: output1 = './output_coordinates.txt'

    if 'SNAPSHOT_MAX' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('SNAPSHOT_MAX')
        snapshot_max = configuration_file['VALUE'][index]

    else: snapshot_max = 95

    if 'SNAPSHOT_MIN' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('SNAPSHOT_MIN')
        snapshot_min = configuration_file['VALUE'][index]

    else: snapshot_min = 2

    if 'CYLINDER_CAT' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('CYLINDER_CAT')
        cylinder_cat = configuration_file['VALUE'][index]

    else:
        cylinder_cat = 'none'
        error += '[+] ERROR: CYLINDER_CAT is a required field!\n'

    if 'HUDF_CAT' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('HUDF_CAT')
        hudf_cat = configuration_file['VALUE'][index]

    else:
        hudf_cat = 'none'
        error += '[+] ERROR: HUDF_CAT is a required field!\n'

    if 'SUBHALO_CAT' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('SUBHALO_CAT')
        subhalo_cat = configuration_file['VALUE'][index]

    else:
        subhalo_cat = 'none'
        error += '[+] ERROR: SUBHALO_CAT is a required field!\n'

    if 'X_CENTER' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('X_CENTER')
        x_center = configuration_file['VALUE'][index]

    else:
        x_center = 2048

    if 'Y_CENTER' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('Y_CENTER')
        y_center = configuration_file['VALUE'][index]

    else:
        y_center = 2048

    if 'X_ANGULAR_FOV' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('X_ANGULAR_FOV')
        x_angular_fov = configuration_file['VALUE'][index]

    else:
        x_angular_fov = 0.0008264462

    if 'Y_ANGULAR_FOV' in configuration_file['ARGS']:
        index = list(configuration_file['ARGS']).index('Y_ANGULAR_FOV')
        y_angular_fov = configuration_file['VALUE'][index]

    else:
        y_angular_fov = 0.0008264462

    return output1, cylinder_cat, hudf_cat, subhalo_cat, snapshot_min, snapshot_max, x_center, y_center, x_angular_fov, y_angular_fov, error



# ======================================================================================================================
#
#       This is the so called 'core' of this script, matching the galaxies observed with galaxies inside the simu-
#       lation and printing it's properties from both the simulation catalog and the Source Extractor output files.
#
# ======================================================================================================================
def main():

    # ------------------------------------------------------------------------------------------------------------------
    #       Parses the filename for the output file with all the objects information
    # ------------------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Performs the galaxy matching algorithm for a lightcone.')
    parser.add_argument('-config', required=True, help='Configuration File')
    args = parser.parse_args()

    # ------------------------------------------------------------------------------------------------------------------
    #       Parses the configuration file
    #
    #       Arguments[0]: OUTPUT1
    #       Arguments[1]: CYLINDER_CAT
    #       Arguments[2]: HUDF_CAT
    #       Arguments[3]: SUBHALO_CAT
    #       Arguments[4]: SNAPSHOT_MIN
    #       Arguments[5]: SNAPSHOT_MAX
    #       Arguments[6]: X_CENTER
    #       Arguments[7]: Y_CENTER
    #       Arguments[8]: X_ANGULAR_FOV
    #       Arguments[9]: Y_ANGULAR_FOV
    #       Arguments[10]: Error
    # ------------------------------------------------------------------------------------------------------------------
    arguments = read_configuration(args.config)

    # ------------------------------------------------------------------------------------------------------------------
    #   Check for errors in the configuration file
    # ------------------------------------------------------------------------------------------------------------------
    if arguments[10] == '':

        # Reads the HUDF catalog
        hudf = ascii.read(arguments[2])

        # constant used to expand the cubes side
        little_h = 0.704

        # ------------------------------------------------------------------------------------------------------------------
        #       Opens the argument-passed file and start printing the information into it
        #
        # ------------------------------------------------------------------------------------------------------------------
        f = open(arguments[0], 'w')

        f.write('# SNAPSHOT\tSUBHALO_ID\tREDSHIFT\tMASS\tX_COORDINATE\tY_COORDINATE\tR_BAND\n')
        f.write('#\tSNAPSHOT\tSimulation snapshot identification number\t[unitless]\t\n' + \
                '#\tSUBHALO_ID\tObjects subhalo identification number inside the simulation\t[unitless]\t\n' +\
                '#\tREDSHIFT\tObjects redshift, it is the same for all the objects within the same cube\t[unitless]\t\n' + \
                '#\tMASS\tObjects mass\t[10**10 Solar]\t\n' +\
                '#\tX_COORDINATE\tX coordinate of the object\t[rad]\t\n' +\
                '#\tY_COORDINATE\tZ coordinate of the object\t[rac]\t\n' + \
                '#\tR_BAND\tR Band absolute magnitude\t[unitless]\t\n' +\
                '#\n')

        # ------------------------------------------------------------------------------------------------------------------
        #       Runs through every cylinder in the lightcone.
        #       Default Range: [02, 95]
        # ------------------------------------------------------------------------------------------------------------------
        for i in range(int(arguments[4]), int(arguments[5])):

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            #       Creates the cylinder and catalog filename and opens it.
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if i < 10:
                catalog_name = arguments[3] + 'cylinder_000' + str(hudf['col1'][i]) + \
                               '_subhalo_catalog.txt'

                cylinder_name = arguments[1] + 'new_noimage_000' + \
                                str(hudf['col1'][i]) + '.fits'
            else:
                catalog_name = arguments[3] + 'cylinder_00' + str(hudf['col1'][i]) + \
                               '_subhalo_catalog.txt'

                cylinder_name = arguments[1] + 'new_noimage_00' + \
                                str(hudf['col1'][i]) + '.fits'

            catalog = ascii.read(catalog_name)
            cylinder = fits.open(cylinder_name)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            #       Load the cube properties.
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            cube_redshift = hudf['col3'][i]

            #fov_angular = cylinder[0].header[13]            # Units: radians
            #fov_linear = cylinder[0].header[22]             # Units: kpc

            # The coordinates the lightcone enters the i-th cube
            ingress_coordinates = np.array([hudf['col10'][i], hudf['col11'][i], hudf['col12'][i]])

            camera_position = np.array([[cylinder[0].header[23]], [cylinder[0].header[24]], [cylinder[0].header[25]]])
            camera_direction = np.array([cylinder[0].header[26], cylinder[0].header[27], cylinder[0].header[28]])
            camera_up = np.array([cylinder[0].header[29], cylinder[0].header[30], cylinder[0].header[31]])

            # Creates a camera object with the previous information
            camera = Camera(camera_position, camera_direction, camera_up)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            #       Runs through every galaxy inside the i-th cube (i.e. cylinder section of the i-th cube).
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            for j in range(0, len(catalog['col1'])):

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                #       Converts the galaxy position from the catalog to a cube with comoving size length 106pc.
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                galaxy_world_position = np.array(
                                       [(catalog['col7'][j]/little_h)*(1/(1+cube_redshift))-ingress_coordinates[0]-float(camera_position[0]),
                                        (catalog['col8'][j]/little_h)*(1/(1+cube_redshift))-ingress_coordinates[1]-float(camera_position[1]),
                                        (catalog['col9'][j]/little_h)*(1/(1+cube_redshift))-ingress_coordinates[2]-float(camera_position[2])
                                       ])

                galaxy_camera_cartesian_position = camera.cameraCoordinates(galaxy_world_position)

                # Geometrically adequate linear field of view
                # fov = ((galaxy_camera_cartesian_position[2]*(little_h*(1+cube_redshift))*(8.2508910903777*(10**(-4))))-25.582322266542)*(1/(little_h*(1+cube_redshift)))

                # 2D Plane coordinates
                y1 = math.atan(galaxy_camera_cartesian_position[0]/galaxy_camera_cartesian_position[2])/(0.5*arguments[8])
                y2 = math.atan(galaxy_camera_cartesian_position[1]/galaxy_camera_cartesian_position[2])/(0.5*arguments[9])

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                #       Cuts the data to consider only the points inside the frame ( within [-1, 1] interval).
                #
                #       --> Use this section if you want to extract more information from the catalog   <-- use j
                #       --> Use this section if you want to extract more information from the hudf file <-- use i
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                if y1 >= -1 and y1 <= 1 and y2 >= -1 and y2 <= 1:
                    f.write(str(catalog['col1'][j]) + '\t' + \
                            str(catalog['col2'][j]) + '\t' + \
                            str(hudf['col3'][i]) + '\t' + \
                            str(catalog['col10'][j]) + '\t' + \
                            str( (y1+1)*arguments[6] ) + '\t' + \
                            str( (y2+1)*arguments[7] ) + '\t' + \
                            str(catalog['col23'][j]) + '\t\n'
                            )
            f.close

    # ------------------------------------------------------------------------------------------------------------------
    #   Exit and print errors
    # ------------------------------------------------------------------------------------------------------------------
    else:
        print('[-] Configuration File Error!\n' + arguments[6])

if __name__ == '__main__':
    main()



