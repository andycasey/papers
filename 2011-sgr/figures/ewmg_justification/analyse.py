import glob, re
from  numpy import *
from scipy import interpolate, ndimage
import matplotlib.pyplot as plt
from pyraf import iraf
import os


def GKItoEPS(GKIfilename):
    '''
    Converts an outputted IRAF GKI file to an EPS file with the same base name and location
    '''

    # Remove all existing first
    if glob.glob('sgi*.eps'): os.system('rm -f sgi*.eps')
    
    EPSfilename = GKIfilename[0:-3] + 'eps'
    
    # Generate the mosaic
    iraf.gflush()
    iraf.gkimosaic(input=GKIfilename, device='eps', nx=1, ny=1, interactive='no', fill='yes', output='', Stdout=1)
    iraf.gflush()
    EPSfiles = glob.glob('sgi*.eps')
    
    i = 0
    while not len(EPSfiles):

        if (i > 15):
            iraf.gkimosaic(input=GKIfilename, device='eps', nx=1, ny=1, interactive='no', fill='yes', output='', Stdout=1)
            iraf.gflush()
            
        EPSfiles = glob.glob('sgi*.eps')
        i += 1

    # File movements
    os.system('mv ' + EPSfiles[0] + ' ' + EPSfilename)
    os.system('rm -f ' + GKIfilename)
    
    return EPSfilename




iraf.noao()
iraf.onedspec()

gauss_val = 2.938/2. # A/pix

synWavelength = loadtxt('LAMBDA_R20.DAT')


obsWavelength = loadtxt('../blue_giant.txt').transpose()[0]

plt.close('all')
spectrumFiles = glob.glob('*.ASC')

data = []
for spectrumFile in spectrumFiles:
    
    string = re.sub('[T|G|V]', ' ', spectrumFile).split()
    Teff = float(string[0])
    
    if 'P' in string[1]:
        logg, feh = string[1].replace('P', ' ').split()
        
        logg = float(logg)/10.
        feh = float(feh)/10.
    
    elif 'M' in string[1]:
        logg, feh = string[1].replace('M', ' ').split()
        
        logg = float(logg)/10.
        feh = float(feh)/-10.
        
    
    # load in the spectra
    
    flux = loadtxt(spectrumFile)
    
    index1, index2 = (synWavelength.searchsorted(obsWavelength[0]) -1), (synWavelength.searchsorted(obsWavelength[-1]) +1)
    synWavelength = synWavelength[index1:index2]
    flux = flux[index1:index2]
    
    
    f = interpolate.interp1d(synWavelength, flux, kind='linear', copy=False)
    new_spectra = f(obsWavelength)
    
    # Convolve the new mapped spectra with a gaussian filter
    
    convolved_spectra = ndimage.gaussian_filter1d(new_spectra, gauss_val)
    
    os.system('rm -f temp.*')
    os.system('cp headers.txt temp.txt')
    
    tempfile = open('temp.txt', 'w')
    for w, i in zip(obsWavelength, convolved_spectra):
        tempfile.write("%4.11f  %3.5f\n"% (w, i,))
        
    tempfile.close()
    
    iraf.onedspec.rspectext(input='temp.txt', output='temp.fits')#, title=spectrumFile), flux='yes', dtype='linear')
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(obsWavelength, convolved_spectra)    
    """
    plotfile = spectrumFile.split('.')[0]+'.gki'
    
    fitprofs = iraf.noao.onedspec.fitprofs(input='temp.fits', reg='5100 5200', positions='locations.pos', fitbackground='yes', fitpositions='all', fitgfwhm='all', fitlfwhm='all', plotfile=plotfile, background='med(5100,5200,1) med(5100,5200,1)', Stdout=1)
    
    
    
    
    EPSfile = GKItoEPS(plotfile)
                    
    # remove the header crap
    fitprofs = fitprofs[3:]
    
    # Assemble the information for the star object
    center, cont, flux, eqw, core, gfwhm, lfwhm = [list()] * 7
    
    # Who loves list comprehension? I do.
    for fitprof in fitprofs: center, cont, flux, eqw, core, gfwhm, lfwhm = [append(obj, val) for obj, val in zip([center, cont, flux, eqw, core, gfwhm, lfwhm], fitprof.split())]

    center, cont, flux, eqw, core, gfwhm, lfwhm = [[float(token.replace('INDEF', '999')) for token in obj.flat] for obj in [center, cont, flux, eqw, core, gfwhm, lfwhm]]
    
    center, cont, flux, core = map(median, [center, cont, flux, core])
    eqw, gfwhm, lfwhm = map(sum, [eqw, gfwhm, lfwhm])


    print Teff, logg, feh
    """
