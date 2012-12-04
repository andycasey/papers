import numpy as np
import matplotlib.pyplot as plt 


v0 = 220 # km/s
vrot = 180 # km/s

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


fig = plt.figure()
ax = fig.add_subplot(111)

# Thick disk stars
ax.scatter(V[thick_disk], pow(pow(U[thick_disk], 2) + pow(W[thick_disk], 2), 0.5), color='r', marker='+', facecolor='none')

# Low-alpha stars
ax.scatter(V[low_alpha], pow(pow(U[low_alpha], 2) + pow(W[low_alpha], 2), 0.5), color='b', marker='o', facecolor='b')

# High-alpha stars
ax.scatter(V[high_alpha], pow(pow(U[high_alpha], 2) + pow(W[high_alpha], 2), 0.5), color='g', marker='o', facecolor='none')

# Plot solar rotation line
ax.plot([-v0, -v0], [0, 400], 'k:', zorder=-1)

# Plot kinematic separation line at Vtot = 180 km/s
x = np.arange(-180, 180)
ax.plot(x, np.sqrt(pow(vrot, 2) - pow(x, 2)), 'k--', zorder=-1)

# Set some limits

ax.set_xlim(-450, 50)
ax.set_ylim(0, 400)
ax.set_xlabel('$V$ (km s$^{-1}$)')
ax.set_ylabel('$(U^{2} + W^{2})^{-1/2}$ (km s$^{-1}$)')

plt.savefig('plot-toombre.eps')
plt.savefig('plot-toombre.pdf')