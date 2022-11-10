# Code to visualize how the 2x2 case evolves per MCMC cycle
import numpy as np
import matplotlib.pyplot as plt

cycles, _, eps, epsSqrd, C_V, __, m, mSqrd, chi = np.loadtxt('../energies.txt', unpack = True)

# plt.plot(cycles, C_V, lw=1, color='red', label='$C_V$')
plt.plot(cycles[50:], chi[50:], lw=1, color='royalblue', label='$\chi$')
# plt.plot(cycles, np.ones(len(cycles)) * .03, ls='dashed', lw=1, color='black', label='$C_V^{exact}$')
plt.plot(cycles[50:], np.ones(len(cycles))[50:] * .008, ls='dashed', lw=1, color='green', label='$\chi^{exact}$')

plt.legend()
plt.show()
