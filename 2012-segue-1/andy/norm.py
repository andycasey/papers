import os
from pyraf import iraf

iraf.noao()
iraf.noao.onedspec()
iraf.obsutil()
iraf.imred()
iraf.ccdr()
iraf.ech()

def normalise_mike_spectra(filename):
    
    os.system('rm -f red.fits red_multic.fits red_multicim.fits red_multicom.fits red_multicimcom.fits endred.fits a.fits b.fits c.fits d.fits e.fits f.fits *.par')
    iraf.noao.onedspec.scopy(input=filename, output="red.fits", aperture="", bands=2, beams="1", format="multispec")
    iraf.noao.onedspec.cont(input="red.fits", output="red_multic.fits", order=6, overrid="yes", markrej="no")
    iraf.images.imutil.imarith(operand1="red.fits", op="/", operand2="red_multic.fits", result="red_multicim.fits")
    iraf.noao.onedspec.scombine("red.fits", "red_multicom.fits", combine="sum", group="images")
    iraf.noao.onedspec.scombine("red_multicim.fits", "red_multicimcom.fits", combine="sum", group="images")
    iraf.images.imutil.imarith(operand1="red_multicom.fits", op="/", operand2="red_multicimcom.fits", result="endred.fits")
    os.system("cp endred.fits %s_norm.fits" % (filename.replace('.fits', ''), ))
    iraf.noao.onedspec.wspectext(input="%s_norm.fits" % (filename.replace('.fits', ''), ), output="%s_norm.txt" % (filename.replace('.fits', ''), ), header="no")
    os.system('rm -f red.fits red_multic.fits red_multicim.fits red_multicom.fits red_multicimcom.fits endred.fits a.fits b.fits c.fits d.fits e.fits f.fits')
    
    return True