#!/bin/bash

# run*.sh
# example LSF submission script for cylinderscript*.pro cylinder extraction controller


module load math/IDL-8.0
bsub -q sunrise -n 4 -R "span[ptile=4] && mem > 80000" -J cylscript -M 88000000 -oo cylscript.out -eo cylscript.err -C 0 idl ~/cosmo_stuff/NatureLetter/BIG_UDF_FILES/cylinderscript_xyz_catsonly
#bsub -q nancy -n 8 -R "span[ptile=8] && mem > 30000" -J cylscript -oo cylscript.out -eo cylscript.err -C 0 idl cylinderscript_75mpc

