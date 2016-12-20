import h5py
import astropy
import astropy.io.ascii as ascii
import numpy as np
import sys
import struct

if __name__=="__main__":
    cat1f = "Illustris-2_RADEC_hudfwide_75Mpc_7_6_yxz_corners.txt"
    cat1d = ascii.read(cat1f)

    types = ['I','L','D','D','D','D','D','D','D','D','D',
             'F','F','F','F','F','F','D','F','I','F','F',
             'F','F','F','F','F','F','F','F','F','F','F',
             'F','F','F','F','F']
    
    print cat1d
    print sys.byteorder
    
    with open('Illustris-2_lightcone_yxz.dat','wb',0) as fOut:
        #print out list of data types/schema for DB table construction
        exrow = cat1d[5]
        for it,exrec in enumerate(exrow.as_void()):
            #print type(exrec), types[it]
            if types[it]=='I':
                print it+1,'32, int'
            elif types[it]=='L':
                print it+1,'64, int'
            elif types[it]=='F':
                print it+1,'32, float'
            elif types[it]=='D':
                print it+1,'64, float'
                    
        for row in cat1d:
            for it,rec in enumerate(row.as_void()):
                if types[it]=='I':
                    fOut.write(np.int32(rec))
                elif types[it]=='L':
                    fOut.write(np.int64(rec))
                elif types[it]=='F':
                    fOut.write(np.float32(rec))
                elif types[it]=='D':
                    fOut.write(np.float64(rec))
                


