import twadynamics as twa
import numpy as np
import matplotlib
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

# %% Plot max of the cavity occupation value as a function of the cavity-vibration coupling strength and the cavity
# frequency

# define cavity vibration system
system = twa.evo()

# define a range of cavity freqency values in cm-1
wavenumberToHartree = 4.55633e-6
freqRange = np.linspace(200, 2400, 200)
freqRangeHa = freqRange * wavenumberToHartree

# define a range of coupling strength values in Hartrees
gRange = np.linspace(0.01e-3, 1e-3, 150)

# True if figure should be saved
savefig = False

# create a array for the heat map of max cavity occupation values
nCavMax = np.zeros((freqRangeHa.size, gRange.size))

# calculate time range to find max value over
totalTime = 16
steps = 8000
timeRange = np.linspace(0, totalTime, steps)
timeRangeAU = timeRange / 2.41889E-5

# set plot colors
colors = plt.cm.plasma(np.linspace(0.5, 1, freqRange.size))

i = 0
for freq in freqRangeHa:
    j = 0
    system.freqC = freq
    for gval in gRange:
        system.g = gval
        cavOc, vibOc, fieldOc, bathOc = system.occ(timeRangeAU, field=False, bath=False)
        nCavMax[i, j] = np.amax(cavOc)
        j += 1
    i += 1

plt.xlabel(r'$g$ / mHa', fontsize=20)
plt.ylabel(r'wavenumber / $\mathrm{cm}^{-1}$', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.imshow(nCavMax, cmap='plasma', aspect='auto',
           extent=[gRange[0] * 1000, gRange[-1] * 1000, freqRange[0], freqRange[-1]], origin='lower', vmax=2)
cbar = plt.colorbar()
cbar.set_label(label=r'max{$\langle N_\mathrm{cav} \rangle$}', size=20)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(16)
if savefig:
    plt.savefig('g_vs_hr.png', dpi=300, bbox_inches='tight')
