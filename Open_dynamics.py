import twadynamics as twa
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('Qt5Agg')

# cavity frequency
omegaC = 1600

# True if figure should be saved
savefig = False

# define on resonant cavity vibration system
system = twa.evo(jccoupling=0.01e-3, huangrhys=0.5, cavityfreq=1590)

# calculate time range
totalTime = 25
steps = 1000
timeRange = np.linspace(0, totalTime, steps)
timeRangeAU = timeRange / 2.41889E-5

qFact = 80
fieldDamp = 1 / (2 * np.pi * qFact)
# Sets the vibrational decay lifetime in ps.
vibLife = 0.2
omegaCHz = 2 * np.pi * 2.998*10**10 * omegaC
qFactBath = omegaCHz * vibLife * 10**-12
bathDamp = 1 / (2 * np.pi * qFactBath)
# bathDamp = 3.92417e-6**2 * (500/0.004) / system.freqV
# bathDamp = 1e-6**2 * (1000/0.004) / system.freqV

system.addbath(bathsize=600, fieldsize=600, bdamp=bathDamp, fdamp=fieldDamp)
# spec = system.spectra(timeRangeAU)
#
# plt.plot(spec)
#
# plt.figure()
nCav, nVib, nField, nBath = system.occ(timeRangeAU)
#
# plt.plot(timeRange, nCav, label='vibration', color='black')
# plt.plot(timeRange, nVib, label='vibration')
# plt.plot(timeRange, nField, label='bath')
# system.addbath(bathsize=500, fieldsize=0, bdamp=3, fdamp=0, dfreq=0.004)
# nCav, nVib, nField, nBath = system.occ(timeRangeAU)

plt.plot(timeRange, nCav, label='cavity', color='black')
plt.plot(timeRange, nVib, label='vibration')
plt.plot(timeRange, nField, label='Field')
plt.plot(timeRange, nBath, label='bath')

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
plt.show()
