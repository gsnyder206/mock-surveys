import math
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import scipy.ndimage
import scipy as sp

#little_f: helper function for arcsinh scaling

def little_f(x, minx, maxx,Q,alph):
	
	f = sp.zeros_like(x)
	
	
	maxx = minx + sp.sinh(Q)/(alph*Q) ; #print maxx
	f = sp.arcsinh(alph*Q*(x-minx))/Q
	f[sp.where(x < minx)]=0.0
	f[sp.where(x > maxx)]=1.0
	
	return f



#strategy for creating nice RGB images, from Lupton et al.:
#  set Q=very small, 1.0e-12
#  adjust alph until you can see the faint features you want and the appropriate level of noise (if desired)
#  keeping alph fixed, adjust Q~1 ish until the bright features are no longer saturated



#Function make_nasa renders images using settings that resemble the outputs of, e.g., NASA and STScI
# press-release images.
#The salient feature is that bright pixels will appear white, giving a natural appearance to galaxies.

#This function both creates the RGB cube and saves the figure.  Most likely you don't want this and should use "make_nasa_interactive"

#inputs:  b,g,r images, arcsinh scaling parameters alpha and Q, and optional noise/PSF parameters
#outputs:  3xNxN array containing scaled RGB values appropriate for passing to matplotlib's imshow function.

def make_nasa(b,g,r,filename,alph,Q,inches=5.0,dpi=72,fwhm_pixels=[0.0,0.0,0.0],sigma_tuple=[0.0,0.0,0.0],zlabel=-1, use_inches=False):


	b = b*1.0
	g = g*1.0
	r = r*1.0

	if fwhm_pixels[0] > 1.0e-5:
		sR=np.zeros_like(r) ; sG = np.zeros_like(g) ; sB = np.zeros_like(b)
		fwhm = fwhm_pixels #pixels, 0.5kpc/pixel
		sigma = fwhm/(2.0*math.sqrt(2.0*math.log(2.0)))
		resR = sp.ndimage.filters.gaussian_filter(r,sigma[0],output=sR)
		resG = sp.ndimage.filters.gaussian_filter(g,sigma[1],output=sG)
		resB = sp.ndimage.filters.gaussian_filter(b,sigma[2],output=sB)

		b = sB
		g = sG
		r = sR

	#the idea is to add sky shot noise *here*, after the sources have been convolved
	if sigma_tuple[0] > 1.0e-8:
		b = b + sigma_tuple[0]*np.random.standard_normal(b.shape)

	if sigma_tuple[1] > 1.0e-8:
		g = g + sigma_tuple[1]*np.random.standard_normal(g.shape)

	if sigma_tuple[2] > 1.0e-8:
		r = r + sigma_tuple[2]*np.random.standard_normal(r.shape)

	b[sp.where(b <= 0.0)]=0.0 ; g[sp.where(g <= 0.0)]=0.0 ; r[sp.where(r <= 0.0)]=0.0


	
	I = (b+g+r)/3.0 + 1.0e-20
	
	
	minval = 0.0
	maxval = np.max(I)
	

	R = little_f(r,minval,maxval,Q,alph)
	G = little_f(g,minval,maxval,Q,alph)
	B = little_f(b,minval,maxval,Q,alph)
	
	imarray = np.asarray([R,G,B])
	
	maxrgbval = np.amax(imarray, axis=0)
	
	changeind = np.where(maxrgbval > 1.0)
	R[changeind] = R[changeind]/maxrgbval[changeind]
	G[changeind] = G[changeind]/maxrgbval[changeind]
	B[changeind] = B[changeind]/maxrgbval[changeind]
	

		
	ind = sp.where(I < 1.0e-10)
	R[ind]=0.0 ; G[ind]=0.0 ; B[ind]=0.0


	imarray = np.asarray(np.transpose([R,G,B]))

	leny = float( len(R[0,:]))
	lenx = float( len(R[:,0]))
	
	inx = lenx/float(dpi)
	iny = leny/float(dpi)
	
	if use_inches==True:
		inx = inches
		iny = inches


	f1 = pyplot.figure(figsize=( inx, iny ), dpi=dpi, frameon=False)
	pyplot.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0, hspace=0.0, wspace=0.0)


	axi=pyplot.axes([0.0,0.0,1.0,1.0], frameon=False)
	axi.set_axis_off()
	axi.imshow(imarray[:,:,:],aspect='auto',interpolation='Nearest')

	if zlabel != -1:
		axi.annotate(str(zlabel),[0.7,0.9])

	
	f1.savefig(filename, dpi=dpi, format='pdf', pad_inches=0)

	pyplot.close(f1)


	
	return imarray




#Function make_general uses the Lupton et al. (2006) approach for rendering CCD images
#The salient feature is that, for 3-filter images, bright pixels will be rendered with a quantitatively correct color in RGB space.

#This function both creates the RGB cube and saves the figure.  Most likely you don't want this and should use "make_interactive"

#inputs:  b,g,r images, arcsinh scaling parameters, and optional noise/PSF parameters
#outputs:  3xNxN array containing scaled RGB values appropriate for passing to matplotlib's imshow function.

def make_general(b,g,r,filename,alph,Q,inches=5.0,dpi=72,fwhm_pixels=0.0,sigma_tuple=[0.0,0.0,0.0],zlabel=-1, use_inches=False):

	fpar = open(filename+'-rgbparams.txt','w')
	fpar.write(filename+'\n')
	fpar.write('alph= {:10e}, Q= {:10e}, inches= {:12.4f}, dpi= {:04d}, fwhm_pixels= {:12.4f}'.format(alph,Q,inches,dpi,fwhm_pixels)+'\n')
	fpar.close()
	
	b = b*1.0
	g = g*1.0
	r = r*1.0

	if fwhm_pixels > 1.0e-5:
		sR=np.zeros_like(r) ; sG = np.zeros_like(g) ; sB = np.zeros_like(b)
		fwhm = fwhm_pixels #pixels, 0.5kpc/pixel
		sigma = fwhm/(2.0*math.sqrt(2.0*math.log(2.0)))
		resR = sp.ndimage.filters.gaussian_filter(r,sigma,output=sR)
		resG = sp.ndimage.filters.gaussian_filter(g,sigma,output=sG)
		resB = sp.ndimage.filters.gaussian_filter(b,sigma,output=sB)

		b = sB
		g = sG
		r = sR

	if sigma_tuple[0] > 1.0e-8:
		b = b + sigma_tuple[0]*np.random.standard_normal(b.shape)

	if sigma_tuple[1] > 1.0e-8:
		g = g + sigma_tuple[1]*np.random.standard_normal(g.shape)

	if sigma_tuple[2] > 1.0e-8:
		r = r + sigma_tuple[2]*np.random.standard_normal(r.shape)

	b[sp.where(b <= 0.0)]=0.0 ; g[sp.where(g <= 0.0)]=0.0 ; r[sp.where(r <= 0.0)]=0.0
	
	I = (b+g+r)/3.0 + 1.0e-20
	
	
	minval = 0.0
	maxval = np.max(I)
	
	
	factor = little_f(I,minval,maxval,Q,alph)/I
	
	R = r*factor
	G = g*factor
	B = b*factor

	
	imarray = np.asarray([R,G,B])
	
	maxrgbval = np.amax(imarray, axis=0)
	
	changeind = np.where(maxrgbval > 1.0)
	R[changeind] = R[changeind]/maxrgbval[changeind]
	G[changeind] = G[changeind]/maxrgbval[changeind]
	B[changeind] = B[changeind]/maxrgbval[changeind]
	
		
	ind = sp.where(I < 1.0e-10)
	R[ind]=0.0 ; G[ind]=0.0 ; B[ind]=0.0


	imarray = np.asarray(np.transpose([R,G,B]))

	leny = float( len(R[0,:]))
	lenx = float( len(R[:,0]))
	inx = lenx/float(dpi)
	iny = leny/float(dpi)

	if use_inches==True:
		inx = inches
		iny = inches


	
	f1 = pyplot.figure(figsize=( inx, iny ), dpi=dpi, frameon=False)
	pyplot.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0, hspace=0.0, wspace=0.0)

	axi=pyplot.axes([0.0,0.0,1.0,1.0], frameon=False)
	axi.set_axis_off()
	axi.imshow(imarray[:,:,:],aspect='auto',interpolation='nearest')

	if zlabel != -1:
		axi.annotate(str(zlabel),[0.7,0.9])

	
	f1.savefig(filename, dpi=dpi, format='png', pad_inches=0)

	pyplot.close(f1)


	
	return imarray








#Function make_interactive uses the Lupton scheme, but only creates the RGB cube for later rendering (i.e., for custom plots)
#inputs:  b,g,r images, arcsinh scaling parameters, and optional noise/PSF parameters
#outputs:  3xNxN array containing scaled RGB values appropriate for passing to matplotlib's imshow function.

def make_interactive(b,g,r,alph,Q,inches=5.0,dpi=72,fwhm_pixels=0.0,sigma_tuple=[0.0,0.0,0.0],zlabel=-1):


	b = b*1.0
	g = g*1.0
	r = r*1.0

	if fwhm_pixels > 1.0e-5:
		sR=np.zeros_like(r) ; sG = np.zeros_like(g) ; sB = np.zeros_like(b)
		fwhm = fwhm_pixels #pixels, 0.5kpc/pixel
		sigma = fwhm/(2.0*math.sqrt(2.0*math.log(2.0)))
		resR = sp.ndimage.filters.gaussian_filter(r,sigma,output=sR)
		resG = sp.ndimage.filters.gaussian_filter(g,sigma,output=sG)
		resB = sp.ndimage.filters.gaussian_filter(b,sigma,output=sB)

		b = sB
		g = sG
		r = sR

	if sigma_tuple[0] > 1.0e-8:
		b = b + sigma_tuple[0]*np.random.standard_normal(b.shape)

	if sigma_tuple[1] > 1.0e-8:
		g = g + sigma_tuple[1]*np.random.standard_normal(g.shape)

	if sigma_tuple[2] > 1.0e-8:
		r = r + sigma_tuple[2]*np.random.standard_normal(r.shape)

	b[sp.where(b <= 0.0)]=0.0 ; g[sp.where(g <= 0.0)]=0.0 ; r[sp.where(r <= 0.0)]=0.0
	
	I = (b+g+r)/3.0 + 1.0e-20
	
	
	minval = 0.0
	maxval = np.max(I)

	
	factor = little_f(I,minval,maxval,Q,alph)/I
	
	R = r*factor
	G = g*factor
	B = b*factor
	
	
	imarray = np.asarray([R,G,B])
	#print imarray.shape
	
	maxrgbval = np.amax(imarray, axis=0)
	#print maxrgbval.shape
	
	changeind = np.where(maxrgbval > 1.0)
	R[changeind] = R[changeind]/maxrgbval[changeind]
	G[changeind] = G[changeind]/maxrgbval[changeind]
	B[changeind] = B[changeind]/maxrgbval[changeind]
	

		
	ind = sp.where(I < 1.0e-10)
	R[ind]=0.0 ; G[ind]=0.0 ; B[ind]=0.0


	imarray = np.asarray(np.transpose([R,G,B]))

	leny = float( len(R[0,:]))
	lenx = float( len(R[:,0]))
	inx = lenx/float(dpi)
	iny = leny/float(dpi)
	
	return imarray


#Function make_interactive_nasa uses the NASA color scheme but only creates the RGB cube, returning it for use in custom plots.
#inputs:  b,g,r images, arcsinh scaling parameters, and optional noise/PSF parameters
#outputs:  3xNxN array containing scaled RGB values appropriate for passing to matplotlib's imshow function.

def make_interactive_nasa(b,g,r,alph,Q,inches=5.0,dpi=72,fwhm_pixels=[0.0,0.0,0.0],sigma_tuple=[0.0,0.0,0.0],zlabel=-1):

	b = b*1.0
	g = g*1.0
	r = r*1.0

	if fwhm_pixels[0] > 1.0e-5:
		sR=np.zeros_like(r) ; sG = np.zeros_like(g) ; sB = np.zeros_like(b)
		fwhm = fwhm_pixels #pixels, 0.5kpc/pixel
		sigma = fwhm/(2.0*math.sqrt(2.0*math.log(2.0)))
		resR = sp.ndimage.filters.gaussian_filter(r,sigma[2],output=sR)
		resG = sp.ndimage.filters.gaussian_filter(g,sigma[1],output=sG)
		resB = sp.ndimage.filters.gaussian_filter(b,sigma[0],output=sB)

		b = sB
		g = sG
		r = sR

	if sigma_tuple[0] > 1.0e-8:
		print("Adding noise to b image: sigma = {:12.6f}".format(sigma_tuple[0]))
		b = b + sigma_tuple[0]*np.random.standard_normal(b.shape)

	if sigma_tuple[1] > 1.0e-8:
		print("Adding noise to g image: sigma = {:12.6f}".format(sigma_tuple[1]))
		g = g + sigma_tuple[1]*np.random.standard_normal(g.shape)

	if sigma_tuple[2] > 1.0e-8:
		print("Adding noise to r image: sigma = {:12.6f}".format(sigma_tuple[2]))
		r = r + sigma_tuple[2]*np.random.standard_normal(r.shape)

	b[sp.where(b <= 0.0)]=0.0 ; g[sp.where(g <= 0.0)]=0.0 ; r[sp.where(r <= 0.0)]=0.0
	
	I = (b+g+r)/3.0 + 1.0e-20
	
	
	minval = 0.0
	maxval = np.max(I)

	
	factor = little_f(I,minval,maxval,Q,alph)/I
	

	R = little_f(r,minval,maxval,Q,alph)
	G = little_f(g,minval,maxval,Q,alph)
	B = little_f(b,minval,maxval,Q,alph)
	
	
	imarray = np.asarray([R,G,B])
	
	maxrgbval = np.amax(imarray, axis=0)
	
	changeind = np.where(maxrgbval > 1.0)
	R[changeind] = R[changeind]/maxrgbval[changeind]
	G[changeind] = G[changeind]/maxrgbval[changeind]
	B[changeind] = B[changeind]/maxrgbval[changeind]
	
	
		
	ind = sp.where(I < 1.0e-10)
	R[ind]=0.0 ; G[ind]=0.0 ; B[ind]=0.0


	imarray = np.asarray(np.transpose([R,G,B]))

	leny = float( len(R[0,:]))
	lenx = float( len(R[:,0]))
	inx = lenx/float(dpi)
	iny = leny/float(dpi)
	
	return imarray




