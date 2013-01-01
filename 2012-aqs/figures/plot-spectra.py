import matplotlib.pyplot as plt
import numpy as np

from specutils import Spectrum1D

accompaning_texts = [
        "C2225316-14437\t4325 / 1.26 / $-$1.26",
        "C2306265-085103\t4170 / 0.85 / $-$1.17",
        "J221821-183424\t4590 / 0.90 / $-$1.61",
        "J223504-152834\t4590 / 2.15 / $-$0.68",
        "J22381-104126\t5140 / 2.96 / $-$1.45"
        ]

spectrum_filenames = ['c2225316-14437_rest.fits', 'c2306265-085103_rest.fits', 'j221821-183424_rest.fits', 'j223504-152834_rest.fits', 'j22381-104126_rest.fits']

order = [1, 0, 3, 2, 4]

# order them
accompaning_texts_ = []
spectrum_filenames_ = []
for n in order[::-1]:
    accompaning_texts_.append(accompaning_texts[n])
    spectrum_filenames_.append(spectrum_filenames[n])

spectrum_filenames = spectrum_filenames_
accompaning_texts = accompaning_texts_

xlims = [6540, 6600]
#xlims = [6150, 6160]

text_x_offset = 2
text_y_offset = 0.10
text_fontsize = 9

label_fontsize = 10


fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95)
ax = fig.add_subplot(111)

for i, (accompaning_text, spectrum_filename) in enumerate(zip(accompaning_texts, spectrum_filenames)):

    spectrum = Spectrum1D.load(spectrum_filename)

    ax.plot(spectrum.disp, spectrum.flux + i, 'k')
    lhs_text, rhs_text = accompaning_text.split('\t')
    ax.text(xlims[0] + text_x_offset, text_y_offset + 1 + i, lhs_text, fontsize=text_fontsize, horizontalalignment='left')
    ax.text(xlims[1] - text_x_offset, text_y_offset + 1 + i, rhs_text, fontsize=text_fontsize, horizontalalignment='right')


ax.set_xlabel('Wavelength, $\lambda$ (${\AA}$)', fontsize=label_fontsize)
ax.set_ylabel('Flux, $F_\lambda$', fontsize=label_fontsize)

ax.get_yticklabels()[0].set_visible(False)
ax.set_xlim(*xlims)
ax.set_ylim(0, i + 1.5)

plt.savefig('spectra-h-alpha.eps')
plt.savefig('spectra-h-alpha.pdf')
