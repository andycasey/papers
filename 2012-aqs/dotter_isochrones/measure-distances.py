import numpy as np
import random
from scipy.interpolate import interp1d
from scipy import ndimage

from galpy.potential import HernquistPotential, LogarithmicHaloPotential, MiyamotoNagaiPotential, NFWPotential, MWPotential
from galpy.potential import plotRotcurve
from galpy.orbit import Orbit


r0 = 8.27   # kpc, Shoenrich 2012
e_r0 = 0.29 # kpc, Shoenrich 2012
v0 = 239.   # km/s, Shoenrich 2012
e_v0 = 9.   # km/s, Shoenrich 2012
G = 4.302e-3# pc/M (km/s)^2, Universe, The
z0 = 0.025  # ??, Default in galpy
e_z0 = 0.0


# Miyamoto-Nagai Disc Potential 
M_disc = 9.3e10 # Solar masses, Helmi et al 2006
a_d = 6.5       # kpc, Helmi et al 2006
b_d = 0.26      # kpc, Helmi et al 2006

# Hernquist Bulge Potential
M_bulge = 3.4e10    # Solar masses, Helmi et al 2006
c_b = 0.7           # kpc, Helmi et al 2006

# Logarithmic Halo Potential
v_halo = 134.   # km/s, Helmi et al 2006
halo_d = 12.    # kpc, Hemli et al 2006


disc_potential  = MiyamotoNagaiPotential(a=a_d/r0, b=b_d/r0)
halo_potential  = LogarithmicHaloPotential(core=halo_d/r0, normalize=0.27)
bulge_potential = HernquistPotential(a=c_b/r0)

# Uncertainties
teff_err = 125
logg_err = 0.2
feh_err = 0.1


stars = [
    {
        'designation':  'J223811-104126',
        'teff': 5190,
        'teff_err': teff_err,
        'logg': 2.93,
        'logg_err': logg_err,
        'feh': -1.63,
        'feh_err': feh_err,

        "v_hel": -235.7, 
        "e_v_hel": 0.1, 
        
        "pmRA": -25.34, 
        "e_pmRA": 2.1, 
        
        "pmDE": -99.53, 
        "e_pmDE": 2.1,

        "ra": 339.54833, 
        "dec": -10.6915, 
        "e_ra": 0.0,
        "e_dec": 0.0,

#        "E(B-V)": 0.12, # From interstellar Na D 1
#        "e_E(B-V)": 0.015, # Half of (0.01, -0.02)
        "E(B-V)": 0.07,
        "e_E(B-V)": 0.01,

        #"E(B-V)": 0.07, # From Schlegel

        # Photometry
        "J": 10.420,
        "e_J": 0.018,

        "K": 9.902,
        "e_K": 0.0172,

        "isochrone": 'J223811-104126_10gyr-1.63+0.4.iso',
        'color': 'b'
    },
    {
        'designation':  'J223504-152834',
        'teff': 4650,
        'teff_err': teff_err,
        'logg': 2.16,
        'logg_err': logg_err,
        'feh': -0.63,
        'feh_err': feh_err,

        "pmRA": 15.9, 
        "e_pmRA": 2.2, 
        "pmDE": -12.76,
        "e_pmDE": 2.2,

        "ra": 338.76875,
        "e_ra": 0.0,

        "dec": -15.47636, 
        "e_dec": 0.0,

        "v_hel": -169.7, 
        "e_v_hel": 0.1, 
        
        'J':    10.363,
        'e_J': 0.025,

        'K':    9.650,
        'e_K': 0.000,
        
        "E(B-V)": 0.04, # Schlegel
        "e_E(B-V)": 0.01, 

        "isochrone": 'J223504-152834_10gyr-0.63+0.4.iso',
        'color': 'r'
    },
    {
        'designation':  'J221821-183424',
        'teff': 4630,
        'teff_err': teff_err,
        'logg': 0.88,
        'logg_err': logg_err,
        'feh': -1.58,
        'feh_err': feh_err,

        "pmRA": -10.6, 
        "e_pmRA": 2.5, 
        "pmDE": -19.27, 
        "e_pmDE": 2.5,

        "ra": 334.58833, 
        "e_ra": 0.0,
        "dec": -18.57453, 
        "e_dec": 0.0,

        "v_hel": -159.5, 
        "e_v_hel": 0.1, 

        'J':    10.34,
        'e_J': 0.021,
        'K':    9.683,
        'e_K': 0.023,
        
        'E(B-V)': 0.0280,
        "e_E(B-V)": 0.01, 

        "isochrone": 'J221821-183434_10gyr-1.58+0.40.iso',
        'color': 'm'
    },
    {
        'designation':  'C230626-085103',
        'teff': 4225,
        'teff_err': teff_err,
        'logg': 0.85,
        'logg_err': logg_err,
        'feh': -1.13,
        'feh_err': feh_err,

        "v_hel": -221.1, 
        "e_v_hel": 0.1, 

        'J':    10.312,
        'e_J': 0.025,
        'K':    9.47,
        'e_K': 0.018,

        "pmRA": -2.49, 
        "e_pmRA": 2.8, 
        "pmDE": -15.26, 
        "e_pmDE": 2.7,

        "ra": 346.61083, 
        "e_ra": 0.0,
        "dec": -8.85133, 
        "e_dec": 0.0,

        'E(B-V)': 0.0461,
        "e_E(B-V)": 0.01, 

        "isochrone": 'C2306265-085103_10Gyr-1.13+0.4.iso',
        'color': 'g'
    },
    {
        'designation':  'C222531-145437',
        'teff': 4365,
        'teff_err': teff_err,
        'logg': 1.25,
        'logg_err': logg_err,
        'feh': -1.22,
        'feh_err': feh_err,

        "v_hel": -156.4, 
        "e_v_hel": 0.1, 

        'J':  10.341,
        'e_J': 0.022,
        'K':   9.572,
        'e_K': 0.024,
        'E(B-V)': 0.0325,
        'e_E(B-V)': 0.01,

        "pmRA": 3.45, 
        "e_pmRA": 2.1, 
        "pmDE": -14.73, 
        "e_pmDE": 2.2,

        "ra": 336.38208,
        "e_ra": 0.0, 
        "dec": -14.911, 
        "e_dec": 0.0,

        "isochrone": 'C2225316-14437_10gyr-1.22+0.4+0.4.iso',
        'color': 'c'
    },
    ]

fig = plt.figure()
ax = fig.add_subplot(111)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)


num_realisations = 10000
c_eep, c_teff, c_logg, c_B, c_V, c_J, c_Ks, = [0, 2, 3, 6, 7, 10, 12]

invisible = []

monte_carlo_data = np.zeros((len(stars), num_realisations, 21))

for j, star in enumerate(stars):

    isochrone_data = np.loadtxt(star['isochrone'])

    # Convert column logTeff to Teff
    isochrone_data[:, c_teff] = pow(10, isochrone_data[:, c_teff])

    # We don't want unevolved stars getting in the mix
    turnoff_index = np.argmax(isochrone_data[:, c_teff])
    rgb = isochrone_data[turnoff_index:]

    distancesJ = []
    distancesK = []

    # X, Y, Z, U, V, W, Vx, Vy, Vz, Vr, Vt, Vphi, E, Lz
    velocities = np.zeros((num_realisations, 21))

    for i in xrange(num_realisations):

        # Randomly select a teff and logg.
        teff_select = random.gauss(star['teff'], star['teff_err'])
        logg_select = random.gauss(star['logg'], star['logg_err'])

        # Find closest points
        distances_to_track = (((teff_select - rgb[:, c_teff])/star['teff_err'])**2 + ((logg_select - rgb[:, c_logg])/star['logg_err'])**2)**0.5

        closest_idx = np.argmin(distances_to_track)
        closest_teff, closest_logg = rgb[closest_idx, c_teff], rgb[closest_idx, c_logg]

        # Closest photometry
        closest_J = rgb[closest_idx, c_J]
        closest_K = rgb[closest_idx, c_Ks]
        closest_B = rgb[closest_idx, c_B]
        closest_V = rgb[closest_idx, c_V]

        # Calculate distance
        extinction = random.gauss(star['E(B-V)'], star['e_E(B-V)'])

        J0 = random.gauss(star['J'], star['e_J']) - extinction * 0.902
        K0 = random.gauss(star['K'], star['e_K']) - extinction * 0.367

        distanceJ = pow(10, (J0 - closest_J + 5)/5.) / 1000. # kpc
        distanceK = pow(10, (K0 - closest_K + 5)/5.) / 1000. # kpc


        distance = distanceJ
        distances_band = 'J'

        # Create splines for the other columns
        
        # Randomly sample from all the other stellar properties
        ra      = random.gauss(star['ra'], star['e_ra'])
        dec     = random.gauss(star['dec'], star['e_dec'])
        pm_ra   = random.gauss(star['pmRA'], star['e_pmRA'])
        pm_dec  = random.gauss(star['pmDE'], star['e_pmDE'])
        velocity= random.gauss(star['v_hel'], star['e_v_hel'])

        # Randomly sample from other galactic properties
        v0_mc   = random.gauss(v0, e_v0)
        r0_mc   = random.gauss(r0, e_r0)
        z0_mc   = random.gauss(z0, e_z0)

        #ra,dec,Dist,UVel,VVel,WVel,RVel,pmRA,pmDE

        orbit = Orbit(vxvv=[ra, dec, distance, pm_ra, pm_dec, velocity], radec=True, uvw=False, lb=False, vo=v0_mc, ro=r0_mc, zo=z0_mc, solarmotion='schoenrich')

        distancesJ.append(distanceJ)
        distancesK.append(distanceK)

        velocity_field = [
                ra,
                dec,
                distance,
                pm_ra,
                pm_dec,
                velocity,
                orbit.x() * r0,
                orbit.y() * r0,
                orbit.z() * r0,
                orbit.U()[0],
                orbit.V()[0],
                orbit.W()[0],
                orbit.vx()[0] * v0_mc,
                orbit.vy()[0] * v0_mc,
                orbit.vz() * v0_mc,
                orbit.vR() * v0_mc,
                orbit.vT() * v0_mc,
                orbit.vphi()[0] * v0_mc,
                orbit.L()[0][2] * v0 * r0, # Lz, 18
                orbit.E(pot=MWPotential) * pow(v0, 2)/2, # Energy
                pow(np.sum(pow(orbit.L()[0][:2], 2)), 0.5) * v0 * r0 # Lperp
            ]
            
        velocities[i, :] = velocity_field


    monte_carlo_data[j, :, :] = velocities

    x_pos, y_pos = velocities[:, 6], velocities[:, 7]

    ax2.scatter(x_pos, y_pos, marker='o', s=3, edgecolor='none', facecolor=star['color'], alpha=0.3)
    invisible.append(ax2.scatter(0, 0, marker='o', s=30, edgecolor='none', facecolor=star['color'], label=star['designation']))

    # Plot this point on the toombre plot
    col_U, col_V, col_W = (9, 10, 11)
    x_vels, y_vels = velocities[:, col_V], pow(pow(velocities[:, col_U], 2) + pow(velocities[:, col_W], 2), 0.5)
    x_delta, y_delta = 1, 1

    x_bins = np.arange(np.floor(np.min(x_vels) - x_delta), np.ceil(np.max(x_vels) + x_delta) + x_delta, x_delta)
    y_bins = np.arange(np.floor(np.min(y_vels) - y_delta), np.ceil(np.max(y_vels) + y_delta) + y_delta, y_delta)

    H, xedges, yedges = np.histogram2d(x_vels, y_vels, bins=(len(x_bins), len(y_bins)), range=None, normed=False, weights=None)

    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]

    ax.scatter(x_vels, y_vels, marker='o', s=3, edgecolor='none', facecolor=star['color'], alpha=0.3)

    invisible.append(ax.scatter(0, 0, marker='o', s=30, edgecolor='none', facecolor=star['color'], label=star['designation']))
    
    nbins = 101
    xi = np.linspace(-1000, 0, nbins + 1)
    yi = np.linspace(-10, 1000, nbins + 1)

    H, xedges, yedges = np.histogram2d(x_vels, y_vels, bins=(xi, yi))

    # Blur it a little
    H2 = ndimage.gaussian_filter(H, sigma=2.0, order=0)

    # Need to calculate sigma contours
    data = H2.flatten()
    levels = data[np.searchsorted(np.cumsum(H), np.array([0.68, 0.95, 0.95]) * np.sum(data))]

    ax.contourf(xi[1:], yi[1:], H2.T, levels=levels, colors=star['color'], alpha=0.5)
    ax.contour(xi[1:], yi[1:], H2.T, levels=levels, colors='k')



    # Calculate +/- for distances
    distances = np.array(distancesJ)

    median = np.median(distances)

    # Calculate positive error
    higher = np.where(distances > median)[0]
    sample = np.zeros((len(higher) * 2))

    sample[len(sample)/2:] = distances[higher]
    sample[:len(sample)/2] = median - (distances[higher] - median)
    positive_err = np.std(sample)

    # Calculate negative error
    lower = np.where(distances < median)[0]
    sample = np.zeros((len(higher) * 2))

    sample[len(sample)/2:] = distances[lower]
    sample[:len(sample)/2] = (median - distances[lower]) + median
    negative_err = np.std(sample)
    print "Distance (%s, %s): %1.2f +/- (%1.2f, %1.2f)" % (star['designation'], distances_band, median, positive_err, -negative_err)




disc_rotation_max = 180

# Plot solar rotation line
ax.plot([-v0, -v0], [0, 10000], 'k:', zorder=-1)

# Plot kinematic separation line at Vtot = 180 km/s
x = np.arange(-disc_rotation_max, disc_rotation_max + 30,)
ax.plot(x, np.sqrt(pow(disc_rotation_max, 2) - pow(x, 2)), 'k--', zorder=-1)

# Plot some Nissen & Schuster comparison data

U, V, W, kind = np.loadtxt('nissen_schuster_2010.data', delimiter='|', skiprows=9, dtype=str, usecols=(-4, -3, -2, -1), unpack=True)

U = np.array(map(float, U))
V = np.array(map(float, V))
W = np.array(map(float, W))


# Thick disk stars
thick_disk = np.where(kind == 'TD')[0]

# High alpha stars
high_alpha = list(np.where(kind == 'high-alpha')[0])
high_alpha.extend(np.where(kind == '(high-alpha)')[0])

high_alpha = np.array(high_alpha)

low_alpha = list(np.where(kind == 'low-alpha')[0])
low_alpha.extend(np.where(kind == '(low-alpha)')[0])

low_alpha = np.array(low_alpha)

# Thick disk stars
ax.scatter(V[thick_disk], pow(pow(U[thick_disk], 2) + pow(W[thick_disk], 2), 0.5), color='k', marker='+', facecolor='none')

# Low-alpha stars
ax.scatter(V[low_alpha], pow(pow(U[low_alpha], 2) + pow(W[low_alpha], 2), 0.5), color='k', marker='s', facecolor='none')

# High-alpha stars
ax.scatter(V[high_alpha], pow(pow(U[high_alpha], 2) + pow(W[high_alpha], 2), 0.5), color='k', marker='d', facecolor='none')


# Limits and labels
ax.set_xlim(-750, 0)
ax.set_ylim(0, 600)
ax.set_xlabel('$V$ (km s$^{-1}$)')
ax.set_ylabel('$(U^2 + W^2)^{1/2}$ (km s$^{-1}$)')

plt.savefig('U-V-W.eps')

ax2.set_xlim(-10, 10)
ax2.set_ylim(-10, 10)

ax.legend(loc=1, prop={'size':11})
ax2.legend(loc=1, prop={'size':11})

ax2.set_xlabel('$X_{\\rm GC}$ (kpc)')
ax2.set_ylabel('$Y_{\\rm GC}$ (kpc)')

[item.set_visible(False) for item in invisible]

# Load Geneva data
geneva_data = np.loadtxt('geneva-survey-cleaned.csv', delimiter=',')
linblad_data = np.zeros((len(geneva_data), 3))

for i, line in enumerate(geneva_data):

    ra, dec, distance, U, V, W, velocity, pm_ra, pm_dec = line

    distance /= 1000.

    orbit = Orbit(vxvv=[ra, dec, distance, pm_ra, pm_dec, velocity], radec=True, uvw=False, lb=False, vo=v0, ro=r0, zo=z0, solarmotion='schoenrich')

    linblad_data[i, :] = [
        orbit.L()[0][2] * v0 * r0, # Lz
        orbit.E(pot=MWPotential) * pow(v0, 2)/2, # Energy
        pow(np.sum(pow(orbit.L()[0][:2], 2)), 0.5) * v0 * r0 # Lperp
        ]


# Draw a Linblad diagram
invisible = []
fig3 = plt.figure(figsize=(8, 6))
fig3.subplots_adjust(left=0.15, bottom=0.10, right=0.95, top=0.95)
#fig3.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.98)

'''
ax3 = fig3.add_subplot(121)

for j, star in enumerate(stars):
    ax3.scatter(monte_carlo_data[j, :, 18], monte_carlo_data[j, :, 20], marker='o', s=3, edgecolor='none', facecolor=star['color'], alpha=0.2)

# Show Geneva data
ax3.scatter(linblad_data[:, 0], linblad_data[:, 2], marker='o', s=3, edgecolor='none', facecolor='k', alpha=0.2)

ax3.set_xlim(-3000, 3000)
ax3.set_ylim(0, 3000)
ax3.set_xlabel('$L_Z$ (kpc km s$^{-1}$)')
ax3.set_ylabel('$L_{\\rm perp}$ (kpc km s$^{-1}$)')

'''

ax4 = fig3.add_subplot(111)

for j, star in enumerate(stars):
    invisible.append(ax4.scatter(0, 0, marker='o', s=30, edgecolor='none', facecolor=star['color'], label=star['designation']))
    ax4.scatter(monte_carlo_data[j, :, 18], monte_carlo_data[j, :, 19], marker='o', s=3, edgecolor='none', facecolor=star['color'], alpha=0.3)

    H, xedges, yedges = np.histogram2d(monte_carlo_data[j, :, 18], monte_carlo_data[j, :, 19], bins=(1000, 1000), range=None, normed=False, weights=None)
    
    nbins = 103
    xi = np.linspace(-3100, 3100, nbins + 1)
    yi = np.linspace(-1.3e5, -0.7e5, nbins + 1)

    H, xedges, yedges = np.histogram2d(monte_carlo_data[j, :, 18], monte_carlo_data[j, :, 19], bins=(xi, yi))

    # Blur it a little
    H2 = ndimage.gaussian_filter(H, sigma=2.0, order=0)

    # Need to calculate sigma contours
    data = H2.flatten()
    levels = data[np.searchsorted(np.cumsum(H), np.array([0.68, 0.95, 0.95]) * np.sum(data))]

    ax4.contourf(xi[1:], yi[1:], H2.T, levels=levels, colors=star['color'], alpha=0.5)
    ax4.contour(xi[1:], yi[1:], H2.T, levels=levels, colors='k')


#ax4.legend(loc=1, prop={'size':11})
[item.set_visible(False) for item in invisible]

# Show Geneva data
ax4.scatter(linblad_data[:, 0], linblad_data[:, 1], marker='o', s=3, edgecolor='none', facecolor='k', alpha=0.3)
ax4.set_xlim(-2800, 2800)
ax4.set_ylim(-1.3e5, -0.7e5)
ax4.set_yticklabels(['$-$%1.1f' % (np.abs(num) * 1e-5, ) for num in ax4.get_yticks()])

ax4.set_xlabel('$L_Z$ (kpc km s$^{-1}$)')
ax4.set_ylabel('Energy (10$^5$ km$^{2}$ s$^{2}$)')

plt.savefig('Lz-E.eps')
plt.savefig('Lz-E.pdf')
plt.savefig('Lz-E.png')