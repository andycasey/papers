
import matplotlib.pyplot as plt
from numpy import loadtxt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


plt.close('all')

fig = plt.figure()

# Giant
giantAxe = fig.add_subplot(211)
wavelength, flux = loadtxt('blue_giant.txt', unpack=True)

giantAxe.plot(wavelength, flux, 'k-')

"""
# FOR TESTING IMAGE ONLY - WE WANT imagebox TO RETURN our .EPS file from IRAF
from matplotlib._png import read_png
fn = get_sample_data("lena.png", asfileobj=False)
arr_lena = read_png(fn)

imagebox = OffsetImage(arr_lena, zoom=0.4)

# END TEST
#xys = (5150., 0.70)

#giantAxeInset = AnnotationBbox(imagebox, xys, xybox=(5300., 0.9), xycoords='data', \
#                  boxcoords='data', pad=0.)
                  
#giantAxe.add_artist(giantAxeInset)

#giantAxeInset = inset_axes(giantAxe, width=1.0, height=3.0, loc=5, bbox_to_anchor=(0.8, 0.5), bbox_transform=giantAxe.transAxes)
#giantAxeInset.set
#giantAxeInset.imshow(arr_lena)
"""

dwarfAxe = fig.add_subplot(212)
wavelength, flux = loadtxt('blue_dwarf.txt', unpack=True)
dwarfAxe.plot(wavelength, flux, 'k-')
dwarfAxe.set_xlabel(r'Wavelength (\AA{})', fontsize=14)
dwarfAxe.set_ylabel('Counts')
giantAxe.set_ylabel('Counts', fontsize=14)
plt.setp(giantAxe.get_xticklabels(), visible=False)

giantAxe.set_xlim(3630., 5650.)
dwarfAxe.set_xlim(3630., 5650.)
plt.draw()
plt.savefig('spectra.eps')

#giantAxeInset = inset_axes(giantAxe, width=5, height=1, loc=5) # Location 5 = on the right, 7 = center right

