import specutils
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure()

fig.subplots_adjust(left=0.10, bottom=0.1, right=0.95, top=0.95, hspace=0.05)

gs = gridspec.GridSpec(2, 1,
                       height_ratios=[1,4]
                       )

diff_ax = plt.subplot(gs[0])
main_ax = plt.subplot(gs[1])

xlims = np.array([4305, 4316])

observed = specutils.Spectrum1D.load('j223811-obs-C.fits')
C_690 = specutils.Spectrum1D.load('j223811-6.90-C.fits')
C_705 = specutils.Spectrum1D.load('j223811-7.05-C.fits')
C_720 = specutils.Spectrum1D.load('j223811-7.20-C.fits')

colors = 'rkb'
C_spectra = [C_690, C_705, C_720]

# Plot main
main_ax.scatter(observed.disp, observed.flux, facecolor='k', s=1)

for color, spectrum in zip(colors, C_spectra):
    main_ax.plot(spectrum.disp, spectrum.flux, c=color, zorder=1 if color != 'k' else 10)

main_ax.plot([0, 5000], [1, 1], ':k')
main_ax.set_xlim(xlims)
main_ax.set_ylim(0, 1.1)

# Plot diff
diff_ax.plot([0, 5000], [0, 0], ':k')

for color, spectrum in zip(colors, C_spectra):

    diff_flux = spectrum.interpolate(observed.disp).flux - observed.flux

    diff_ax.plot(observed.disp, diff_flux, c=color, zorder=1 if color != 'k' else 10)

diff_ax.xaxis.set_visible(False)
diff_ax.set_xlim(xlims)
diff_ax.set_ylim(-0.15, 0.15)
diff_ax.set_yticks([-0.10, 0.0, 0.1])

main_ax.set_xlabel('Wavelength, $\lambda$ (${\AA}$)')
main_ax.set_ylabel('Flux, $F_{\lambda}$')

diff_ax.set_ylabel('$\Delta{F_{\lambda}}$')

main_ax.text(4309, 0.15, 'J223811-104126', color='k')
main_ax.text(4312, 0.15, "[C/Fe] =            0.05", color='k', horizontalalignment='left')
main_ax.text(4312, 0.08, "5190/2.93/$-$1.43", color='k', horizontalalignment='left')
main_ax.text(4313.2, 0.15, '$-$0.10', color=colors[0], horizontalalignment='left')
main_ax.text(4314.85, 0.15, '0.20', color=colors[-1], horizontalalignment='left')

plt.savefig('J223811-C.pdf')
