


####
#### Name:    make_allsubmit_script.py
#### Author:  Greg Snyder   gsnyder@stsci.edu
#### Purpose:  generates LSF submission script for a bunch of Sunrise calculations (Jonsson 2006) in a single directory.
####  

import struct
import numpy as np
import asciitable
import pyfits
import subprocess
import glob
import os
import sys

if __name__=="__main__":
    if (len(sys.argv)!=4):
        print "Error: Requires two strings for file template and one string argument (last is queue)."
        quit()

    scriptfile = "allsubmitscript.sh"
    print 'Making '+scriptfile

    queue = sys.argv[3]
    #int1 = int(sys.argv[3]) ; int2 = int(sys.argv[4])

    f = open(scriptfile,'w')
    f.write('#!/bin/bash\n')

    fileglob = np.sort(np.array(glob.glob(sys.argv[1] + '*'+sys.argv[2])))

    print queue
    print fileglob

    #grid = np.arange(int1, int2+1)
    for file in fileglob:
        #numstring = '{:04d}'.format(num)
        command = 'bsub -q '+queue+' -g /GFS < '+file
        print command ; f.write(command+'\n')

    f.close()
