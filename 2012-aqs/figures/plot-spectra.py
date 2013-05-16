import matplotlib.pyplot as plt
import numpy as np

from specutils import Spectrum1D

accompaning_texts = [
        "C2225316-14437\t\t\t4365 / 1.25 / $-$1.22",
        "C2306265-085103\t\t\t4225 / 0.85 / $-$1.13",
        "J221821-183424\t\t\t4630 / 0.88 / $-$1.58",
        "J223504-152834\t\t\t4650 / 2.16 / $-$0.63",
        "J223811-104126\t\t\t5190 / 2.93 / $-$1.63"
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
xlims = [4835, 4885]

text_x_offset = 2
text_y_offset = 0.10
text_fontsize = 11

label_fontsize = 10


fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95)
ax = fig.add_subplot(111)

for i, (accompaning_text, spectrum_filename) in enumerate(zip(accompaning_texts, spectrum_filenames)):

    spectrum = Spectrum1D.load(spectrum_filename)

    ax.plot(spectrum.disp, spectrum.flux + i, 'k')
    #lhs_text, rhs_text = accompaning_text.split('\t')
    #ax.text(xlims[0] + text_x_offset, text_y_offset + 1 + i, lhs_text, fontsize=text_fontsize, horizontalalignment='left')
    #ax.text(xlims[1] - text_x_offset, text_y_offset + 1 + i, rhs_text, fontsize=text_fontsize, horizontalalignment='right')
    ax.text(xlims[0] + text_x_offset, text_y_offset + i + 1, accompaning_text, fontsize=text_fontsize, horizontalalignment='left')


ax.set_xlabel('Wavelength, $\lambda$ (${\AA}$)', fontsize=label_fontsize)
ax.set_ylabel('Flux, $F_\lambda$', fontsize=label_fontsize)

ax.get_yticklabels()[0].set_visible(False)
ax.set_xlim(*xlims)
ax.set_ylim(0, i + 1.5)

plt.savefig('spectra-sample.eps')
plt.savefig('spectra-sample.pdf')
