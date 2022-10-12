import twadynamics as twa
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('Qt5Agg')

# cavity frequency
omegaC = 1400

# True if figure should be saved
savefig = False

# define on resonant cavity vibration system
system = twa.evo(huangrhys=0.5, cavityfreq=omegaC)

# calculate time range
time = np.array([18])

qFact = 6
fieldDamp = 1 / (2 * np.pi * qFact)
qFactBath = 602.7672470628711  # gives vibrational lifetime of 2 ps
bathDamp = 1 / (2 * np.pi * qFactBath)
# bathDamp = 3.92417e-6**2 * (500/0.004) / system.freqV
# bathDamp = 1e-6**2 * (1000/0.004) / system.freqV

system.addbath(bathsize=500, fieldsize=500, bdamp=bathDamp, fdamp=fieldDamp)

# heat map variables:
# Quality factor values
QFRange = np.linspace(2, 80, 10)
# Converts quality factor into the cavity external field coupling factor kappa
fDampRange = 1/(2 * np.pi * QFRange**2)
# cavity-vibration coupling values
gRange = np.linspace(1e-8, 0.1e-03, 10)

percEmiss = np.zeros((gRange.size, QFRange.size))

t = 0
i = 0
for fDamp in fDampRange:
    j = 0
    for g in gRange:
        system.g = g
        system.updatefCoupling(fDamp)
        percEmiss[j, i] = system.occ(time / 2.41889E-5, cav=False, vib=False, bath=False)[2]/0.5 * 100
        t += 1
        print(t)
        j += 1
    i += 1

plt.figure()
plt.ylabel(r'$g$ / mHa', fontsize=20)
plt.xlabel(r'$\mathrm{QF}^{1/2}$', fontsize=20)
# plt.xticks(np.arange(0, max(kappa_range * 10**6)+5, 5))
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.imshow(percEmiss, extent=[QFRange[0], QFRange[-1], gRange[0] * 1e3, gRange[-1] * 1e3], aspect='auto', origin='lower')
cbar = plt.colorbar()
cbar.set_label(label='% emission', size=20)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(16)
if savefig:
    plt.savefig('percent_emission_Qsqr_heat_high_res', dpi=300, bbox_inches='tight')
