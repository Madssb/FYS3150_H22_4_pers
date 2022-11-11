# Code to visualize how many MCMC cycles is needed for good convergence against
# analytical results

import numpy as np
import matplotlib.pyplot as plt

# Analytical values
anaEps = -1.995982
anaM = 0.99866
anaC_V = 0.03208233
anaX = 0.0080282

nCycles, eps, magn, C_V, X = np.loadtxt('../energies.txt', unpack=True)

plt.plot(nCycles, C_V, lw=1)
plt.plot(nCycles, np.ones(len(nCycles)) * anaC_V, lw=1, ls='dashed')
plt.show()
