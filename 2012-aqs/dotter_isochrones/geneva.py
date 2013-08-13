import numpy as np

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

#ra,dec,Dist,UVel,VVel,WVel,RVel,pmRA,pmDE


disc_potential  = MiyamotoNagaiPotential(a=a_d/r0, b=b_d/r0)
halo_potential  = LogarithmicHaloPotential(core=halo_d/r0, normalize=0.27)
bulge_potential = HernquistPotential(a=c_b/r0)


timesteps = np.linspace(0, 100, 1000)
c_ra, c_dec, c_dist, c_U, c_V, c_W, c_rvel, c_pma, c_pmde = xrange(9)

data = np.loadtxt('geneva-survey-cleaned.csv', delimiter=',')

linblad_data = np.zeros((len(data), 6))
measured_uvw = np.zeros((len(data), 3))

for i, line in enumerate(data):

    print i
    ra, dec, distance, U, V, W, velocity, pm_ra, pm_dec = line

    distance /= 1000.

    orbit = Orbit(vxvv=[ra, dec, distance, pm_ra, pm_dec, velocity], radec=True, uvw=False, lb=False, vo=v0, ro=r0, zo=z0, solarmotion='schoenrich')

    measured_uvw[i, :] = [orbit.U()[0], orbit.V()[0], orbit.W()[0]]

    #orbit.integrate(timesteps, pot=[halo_potential, bulge_potential, disc_potential])

    linblad_data[i, :3] = [
        orbit.L()[0][2] * v0 * r0, # Lz
        orbit.E(pot=MWPotential) * pow(v0, 2)/2, # Energy
        pow(np.sum(pow(orbit.L()[0][:2], 2)), 0.5) * v0 * r0 # Lperp
        ]
        
    #orbit.E(),
        


fig = plt.figure()
ax = fig.add_subplot(311)
ax.scatter(data[:, c_U], measured_uvw[:, 0])

limits = ax.get_xlim()
ax.plot(limits, limits, 'k:', zorder=-1)
ax.set_xlim(limits)
ax.set_ylim(limits)
ax.set_xlabel('$U$ (Holmberg)')
ax.set_ylabel('$U$ (Mine)')

ax = fig.add_subplot(312)
ax.scatter(data[:, c_V], measured_uvw[:, 1])
limits = ax.get_xlim()
ax.plot(limits, limits, 'k:', zorder=-1)
ax.set_xlim(limits)
ax.set_ylim(limits)
ax.set_xlabel('$V$ (Holmberg)')
ax.set_ylabel('$V$ (Mine)')



ax = fig.add_subplot(313)
ax.scatter(data[:, c_W], measured_uvw[:, 2])

limits = ax.get_xlim()
ax.plot(limits, limits, 'k:', zorder=-1)
ax.set_xlim(limits)
ax.set_ylim(limits)
ax.set_xlabel('$W$ (Holmberg)')
ax.set_ylabel('$W$ (Mine)')

plt.savefig('geneva-uvw-check.png')


