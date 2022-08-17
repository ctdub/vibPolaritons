import twadynamics as twa
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

# %% Plot a range of cavity frequencies

# define a range of cavity freqency values
wavenumberToHartree = 4.55633e-6
freqRange = np.linspace(1000, 2400, 6)
freqRangeHa = freqRange * wavenumberToHartree

# True if figure should be saved
savefig = False

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
    if i == 0:
        axs[freqRangeHa.size - 1 - i].set_ylim(0, 1.6)
    else:
        axs[freqRangeHa.size - 1 - i].set_ylim(-0.1, 1.6)
    axs[freqRangeHa.size - 1 - i].set_aspect(0.05)
    axs[freqRangeHa.size - 1 - i].tick_params(axis='y', which='major', labelsize=16)
    if i > 0:
        axs[freqRangeHa.size - 1 - i].spines['bottom'].set_visible(False)
        axs[freqRangeHa.size - 1 - i].tick_params(bottom=False)
        axs[freqRangeHa.size - 1 - i].xaxis.set_ticklabels([])
    i += 1

plt.xlabel('time / ps', fontsize=20)
plt.xticks(fontsize=16)
if savefig:
    plt.savefig('ncav-traces-freqdependence-gval=6e-4.png', dpi=300, bbox_inches='tight')

# %% Plot multiple dynamics at different frequencies on top of each other.

# True if figure should be saved
savefig = False

# define a range of cavity freqency values
wavenumberToHartree = 4.55633e-6
freqRange = np.array([400])
freqRangeHa = freqRange * wavenumberToHartree

# define on resonant cavity vibration system with a coupling strength and Huang-Rhys factor
system = twa.evo(1600, 1600, 0.75, 0.6e-3)

# calculate time range
totalTime = 8
steps = 16000
timeRange = np.linspace(0, totalTime, steps)
timeRangeAU = timeRange / 2.41889E-5

# set plot colors
colors = plt.cm.jet(np.linspace(0.5, 1, freqRange.size))

# calculate the time evolution of the cavity and vibration occupation values and plot them
plt.figure(2)
i = 0
for freq in freqRangeHa:
    system.freqC = freq
    cavOc, vibOc = system.isolated(timeRangeAU)
    plt.plot(timeRange, cavOc, color=colors[i], label=freqRange[i])
    i += 1

plt.gca().spines['left'].set_linewidth(2)
plt.gca().spines['bottom'].set_linewidth(2)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('time / ps', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel(r'$\langle N_{x} \rangle$', fontsize=20)
plt.xticks(fontsize=16)
plt.gca().xaxis.set_tick_params(width=2)
plt.gca().yaxis.set_tick_params(width=2)
legend = plt.legend(title='$x = $', loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, fontsize=16)
legend.get_title().set_fontsize('20')
legend.get_frame().set_edgecolor('black')
if savefig:
    plt.savefig('/plots/dynamics_twa_q3.png', dpi=300, bbox_inches='tight')

# %% Plot max of the cavity occupation value as a function of the cavity-vibration coupling strength and the cavity
# frequency

# define a range of cavity freqency values in cm-1
wavenumberToHartree = 4.55633e-6
freqRange = np.linspace(800, 2400, 150)
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
        cavOc, vibOc = system.isolated(timeRangeAU)
        nCavMax[i, j] = np.amax(cavOc)
        j += 1
    i += 1


plt.xlabel(r'$g$ / mHa', fontsize=20)
plt.ylabel(r'wavenumber / $\mathrm{cm}^{-1}$', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.imshow(nCavMax, cmap='plasma', aspect='auto',
           extent=[gRange[0] * 1000, gRange[-1] * 1000, freqRange[0], freqRange[-1]], origin='lower')
cbar = plt.colorbar()
cbar.set_label(label=r'max{$\langle N_\mathrm{cav} \rangle$}', size=20)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(16)
if savefig:
    plt.savefig('g_vs_hr.png', dpi=300, bbox_inches='tight')

#%% plot alone

# True if figure should be saved
savefig = True

plt.xlabel(r'$g$ / mHa', fontsize=20)
plt.ylabel(r'wavenumber / $\mathrm{cm}^{-1}$', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.imshow(nCavMax, cmap='plasma', aspect='auto',
           extent=[gRange[0] * 1000, gRange[-1] * 1000, freqRange[0], freqRange[-1]], origin='lower')
cbar = plt.colorbar()
cbar.set_label(label=r'max{$\langle N_\mathrm{cav} \rangle$}', size=20)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(16)
if savefig:
    plt.savefig('g-vs-cavfreq-heatmap.png', dpi=300, bbox_inches='tight')