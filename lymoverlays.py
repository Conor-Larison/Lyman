import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

from astropy.utils.data import get_pkg_data_filename
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.io import ascii
import sys
from astropy import wcs
import math as mt
import numpy as np
import glob
from scipy.ndimage import filters
from scipy import stats

def main():

	writedir = '/home/conor/Lyman/Project/'
	writedir1 = '/home/conor/Lyman/LymProject/'	


	fits_list = sorted(glob.glob(writedir+'*.fits'))
	fits_list1 = sorted(glob.glob(writedir1+'*.fits'))
	i = 0
	j = 0
	while j<len(fits_list):
		name = fits_list[j].split("_")[0][-6:]


		if name not in fits_list1[i]:
			i += 1
			continue


		fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols = 3)
		
		a = fits.open(fits_list[j])
		image_file = get_pkg_data_filename(fits_list[j])
		fits.info(image_file)

		b = fits.open(fits_list1[i])
		image_file1 = get_pkg_data_filename(fits_list1[i])

		image_data = fits.getdata(image_file)[39:59,39:59]
		image_data1 = fits.getdata(image_file1)[39:59,39:59]
		simg = filters.gaussian_filter(image_data1, .5)
		std = 1.4826 * stats.median_absolute_deviation(simg)
		std = np.median(std)

		ax1.imshow(image_data, cmap='binary',origin='lower')
		m1 = ax1.imshow(image_data, cmap='binary',origin='lower')
		cbar = fig.colorbar(m1, ax=ax1,shrink = .3)
		cbar.ax.tick_params(labelsize=7.5) 
		min = np.amin(simg)
		max = np.amax(simg)
		ax1.set_title('IR Stamp',fontsize = 8)


		ax2.imshow(image_data1, cmap='binary',origin='lower')
		m2 = ax2.imshow(image_data1, cmap='binary',vmin = min,vmax = max,origin='lower')
		cbar2 = fig.colorbar(m2,ax=ax2,shrink = .3)
		cbar2.ax.tick_params(labelsize=7.5)
		ax2.set_title('Lyman-alpha Stamp',fontsize = 8)
		

		ax3.imshow(image_data, cmap='binary',origin='lower')
		m3 = ax3.imshow(image_data, cmap='binary',origin='lower')
		cbar3 = fig.colorbar(m3,ax=ax3,shrink = .3)
		cbar3.ax.tick_params(labelsize=7.5)
		ax3.contour(simg,[std,2*std,3*std,4*std,5*std,6*std], colors='red',alpha =.3)
		ax3.set_title('IR With Lyman Contours', fontsize = 8)

		ax1.tick_params(labelsize=8.5)
		ax2.tick_params(labelsize=8.5)
		ax3.tick_params(labelsize=8.5)


		j+=1
		
		pp = PdfPages('/home/conor/Lyman/LymCatalog/model_'+name+'.pdf')
		pp.savefig(fig)
		pp.close()



main()