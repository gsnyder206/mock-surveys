import os
import glob
import numpy as np
import string


if __name__=="__main__":
    files = np.asarray(glob.glob('*.fits'))

    
    for fn in files:

        fieldstr=fn[-36:-18]
        endstr = '_v1_lightcone.fits'

        oldend=fieldstr+endstr
        
        newfn=fn.rstrip(oldend).lower()+'_'+(fieldstr.lower()).replace('_','-')+endstr
        
        print(fn)
        print(newfn)
        print('')

        os.rename(fn,newfn)
        
