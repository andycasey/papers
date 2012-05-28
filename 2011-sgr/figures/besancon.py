import matplotlib.pyplot as plt
import glob
import atpy
from numpy import *
import re
import ephem

filenames = glob.glob('*.resu')

bandwidth = 10.0
x = arange(-400., 400., 1.)


besanconVgsr = []
color = []
mag = []
VgsrCorrection = 0.0
for filename in filenames:
    
    datafile = open(filename, 'r')
    content = datafile.readlines()
    
    for line in content:
        if re.search('a\(2000\.0\)\s+:\s+', line):
            
            # could be done better, but fuck it
            a, b, c = line.split(':')
            b, a = b.split('Step')
            
            ra1, ra2 = b.split('-')
            ra1 = ra1.replace('h', ':').replace('m', ':').strip(' s').replace(' ', '')
            ra2 = ra2.replace('h', ':').replace('m', ':').strip(' s').replace(' ', '')
            ra1 = ephem.Equatorial(ra1, 0., epoch=ephem.J2000)
            ra2 = ephem.Equatorial(ra2, 0., epoch=ephem.J2000)
            
            ra = 0.5*(ra2.ra - ra1.ra) + ra1.ra
            
            c, d = c.split('Step')
            dec1, dec2 = c.split(' - ')
            
            dec1 = dec1.replace('d', ':').replace('m', ':').strip(' s').replace(' ', '')
            dec2 = dec2.replace('d', ':').replace('m', ':').strip(' s').replace(' ', '')
            dec1 = ephem.Equatorial(0., dec1, epoch=ephem.J2000)
            dec2 = ephem.Equatorial(0., dec2, epoch=ephem.J2000)
            
            dec = 0.5*(dec2.dec - dec1.dec) + dec1.dec
            
            coords = ephem.Equatorial(ra, dec, epoch=ephem.J2000)
            galactic = ephem.Galactic(coords, epoch=1958)
            l, b = [degrees(_) for _ in [galactic.long, galactic.lat]]

            # Get heliocentric and galactocentric corrections (Website stefan gave me)
            #VgsrmVobs = 220.0*sin(l)*cos(b) + 9.0*cos(l)*cos(b) + 12.0*sin(l)*cos(b) + 7.0*sin(b) 
      
            # Yanny 2009
            #VgsrmVobs = 10.1*cos(b)*cos(l) + 224.0*cos(b)*sin(l) + 6.7*sin(b)
    
            # Beaulieu & Freeman, et al 2000 (they call it Vlos,GC)
            VgsrCorrection = 220.0*sin(radians(l))*cos(radians(b)) + 16.5*(sin(radians(b))*sin(radians(25)) + cos(radians(b))*cos(radians(25))*cos(radians(l-53)))
 
            
        elif re.search('\(l =\s+\d{1,3}\.\d; b =\s+\d{1,2}\.\d\)', line):
            
            l, b, d = line.split(';')
            
            l = float(l.strip(' =l('))
            b = float(b.split(')')[0].strip(' b='))
                      
            VgsrCorrection = 220.0*sin(radians(l))*cos(radians(b)) + 16.5*(sin(radians(b))*sin(radians(25)) + cos(radians(b))*cos(radians(25))*cos(radians(l-53)))

            
        elif re.search('(-?\d{1,3}\.?\d{0,5}\s+){24}', line):
            besanconVgsr.append(float(line.split()[15]) + VgsrCorrection)
            color.append(float(line.split()[8])) #B-V U-B V-I V-K V
            mag.append(float(line.split()[1])) # Mv
            

y = zeros(len(x))
for vgsr in besanconVgsr:
    y += exp((-(x-vgsr)**2)/(2*bandwidth**2))/sqrt(2*pi)
    
 


observations = atpy.Table('combined.xml', verbose=False)


subsetLogic = '(observations.K_Giants == True) & (-250. < observations.Vgsr) & (observations.Vgsr < 250.)'


# Generate the subset

observationsSubset = observations.where(eval(subsetLogic)).data



z = zeros(len(x))
for vgsr in observationsSubset['Vgsr']:
    z += exp((-(x-vgsr)**2)/(2*bandwidth**2))/sqrt(2*pi)

# Scale the kernel
y = y / (len(besanconVgsr) * bandwidth)
z = z / (len(observationsSubset['Vgsr']) * bandwidth)
    

plt.close('all')
fig = plt.figure()

labelsize = 14.5
ax = fig.add_subplot(111)

ax.plot(x, y, 'k--', label=r'Besan\c{c}on') # besancon
ax.plot(x, z, 'k-', label='Observed K-giants') # observed

ax.set_xlabel('$V_{GSR}$ (km s$^{-1}$)', labelsize=labelsize)
ax.set_ylabel('$\sum{}P(V_{GSR})$',labelsize=labelsize)
ax.legend()
plt.draw()
plt.savefig('besancon.eps')

plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(mag, color)
ax.set_xlabel('B-V')
ax.set_ylabel('$M_V$')

plt.draw()
plt.savefig('besancon-cmd.eps')
plt.close('all')
