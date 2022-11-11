# Code to visualize how estimates of <e> and <|m|> evolve with the number of
# MCMC cycles

import numpy as np
import matplotlib.pyplot as plt

nCycles, eps, magn, C_V, X = np.loadtxt('../unordered_2by2_lattice.txt', unpack=True)

plt.plot(nCycles, eps)
plt.show()
