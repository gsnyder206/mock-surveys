import astropy
import astropy.io.fits as fits
import webbpsf
from webbpsf import wfirst
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np

wfi=wfirst.WFI()

fs=wfi.filter_list

print(fs)
print(wfi.detector)
print(wfi.detector_position)


fig=pyplot.figure(figsize=(8.0,4.0),dpi=300)

fig.subplots_adjust(left=0.0,bottom=0.0,top=1.0,right=1.0,hspace=0.0,wspace=0.0)

for i,fname in enumerate(fs):
    wfi.filter=fname
    this_psf=wfi.calc_psf()
    this_psf.writeto('WFI_'+fname+'_SCA01_center_psf_phaseb.fits',overwrite=True)
    ax=fig.add_subplot(2,4,i+1)
    ax.set_axis_off()
    ax.imshow(np.log10(this_psf[0].data))
    ax.annotate(fname,(0.5,0.15),xycoords='axes fraction',ha='center',va='center',color='white',size=16)

fig.savefig('WFI_PSFs.pdf',dpi=300)
