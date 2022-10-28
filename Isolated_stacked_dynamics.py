import twadynamics as twa
import numpy as np
import matplotlib
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

# %% Plot a range of cavity frequencies

# define a range of cavity freqency values
wavenumberToHartree = 4.55633e-6
freqRange = np.linspace(450, 2400, 6)
freqRangeHa = freqRange * wavenumberToHartree

# True if figure should be saved
savefig = True

# define on resonant cavity vibration system
system = twa.evo(jccoupling=0.6e-3)

# calculate time range
totalTime = 0.4
steps = 1000
timeRange = np.linspace(0, totalTime, steps)
timeRangeAU = timeRange / 2.41889E-5

# set plot colors
colors = plt.cm.jet(np.linspace(0.5, 1, freqRange.size))

# calculate the time evolution of the cavity and vibration occupation values and plot them on top of each other
fig, axs = plt.subplots(freqRangeHa.size)
i = 0
for freq in freqRangeHa:
    system.freqC = freq
    cavOc, vibOc, fieldOc, bathOc = system.occ(timeRangeAU)
    axs[freqRangeHa.size - 1 - i].plot(timeRange, cavOc, color=colors[i])
    axs[freqRangeHa.size - 1 - i].spines['right'].set_visible(False)
    axs[freqRangeHa.size - 1 - i].spines['top'].set_visible(False)
    axs[freqRangeHa.size - 1 - i].spines['left'].set_visible(False)
    axs[freqRangeHa.size - 1 - i].set_yticks([0, system.hr])
    axs[freqRangeHa.size - 1 - i].set_yticklabels(['0', '0.75'])
    axs[freqRangeHa.size - 1 - i].set_xticks(np.round(np.linspace(0, totalTime, 5, endpoint=True), 2))
    if i == 0:
        axs[freqRangeHa.size - 1 - i].set_ylim(0, 2)
    else:
        axs[freqRangeHa.size - 1 - i].set_ylim(-0.1, 2)
    axs[freqRangeHa.size - 1 - i].set_aspect(0.04)
    axs[freqRangeHa.size - 1 - i].tick_params(axis='y', which='major', labelsize=16)
    if i > 0:
        axs[freqRangeHa.size - 1 - i].spines['bottom'].set_visible(False)
        axs[freqRangeHa.size - 1 - i].tick_params(bottom=False)
        axs[freqRangeHa.size - 1 - i].xaxis.set_ticklabels([])
    i += 1



plt.xlabel('time / ps', fontsize=20)
plt.xticks(fontsize=16)
if savefig:
    plt.savefig('ncav-traces-freqdependence-gval=6e-4_450-2400freq.png', dpi=300, bbox_inches='tight')

