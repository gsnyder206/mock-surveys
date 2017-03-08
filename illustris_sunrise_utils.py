import os
import sys
import numpy as np
import glob


#this function takes as inputs the return values of illustris_api_utils.get_subhalo()
#this is a snapshot file and subhalo metadata object

def setup_sunrise_illustris_subhalo(snap_cutout,subhalo_object,verbose=True,clobber=True):

    snap_dir = os.path.dirname(os.path.abspath(snap_cutout))
    
    print("Setting up sunrise run in.. ", snap_dir)


    
    return

