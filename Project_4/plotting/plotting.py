# Code to visualize how estimates of <e> and <|m|> evolve with the number of
# MCMC cycles

import numpy as np
import matplotlib.pyplot as plt

def analyticalValues():
    '''
    Return the analytical values of energy and magnetization per spin,
    specific heat cap. and susceptibility for 2x2 case
    '''

    # Expected values for <e>, <e^2>, <|m|> and <m^2>
    e = -2 * np.sinh(8) / (np.cosh(8) + 3)
    e_sqrd = 4 * np.cosh(8) / (np.cosh(8) + 3)
    m = 1. / 2 * (np.exp(8) + 2) / (np.cosh(8) +3)
    m_sqrd = 4 * (np.exp(8) + 1) / (np.cosh(8) + 3)

    return e, e_sqrd, m, m_sqrd

def plotConvergence(filename, L, include_analytical=None):
    '''
    Take a .txt file with computed values and plot the
    evolment with number of MCMC cycles
    '''

    nCycles, eps, magn, C_V, X = np.loadtxt(filename, unpack=True)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10.6, 6), sharex=True)

    fig.suptitle(f'{L}x{L} lattice')

    ax1.set_title('Energy per spin')
    ax1.plot(nCycles, eps, lw=1, color='red', label='$\\overline{\\varepsilon}$')
    ax1.set_ylabel('[J]')

    ax2.set_title('Magnetization per spin')
    ax2.plot(nCycles, magn, lw=1, color='royalblue', label='$\\overline{m}$')
    ax2.set_ylabel('M')

    ax3.set_title('Specific heat capacity')
    ax3.plot(nCycles, C_V, lw=1, color='black', label='$\\overline{C_V}$')
    ax3.set_ylabel('$C_V$')
    ax3.set_xlabel('No. of cycles')

    ax4.set_title('Susceptibility')
    ax4.plot(nCycles, X, lw=1, color='green', label='$\\overline{\chi}$')
    ax4.set_ylabel('$\chi$')
    ax4.set_xlabel('No. of cycles')

    if(include_analytical == True):

        e, e_sqrd, m, m_sqrd = analyticalValues()

        y = np.ones(len(nCycles))

        ax1.plot(nCycles, y * e, color='green', lw=1, ls='dashed', label='$<\\varepsilon>$')
        ax2.plot(nCycles, y * m, color='black', lw=1, ls='dashed', label='<m>')

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()

    plt.tight_layout()


plotConvergence('../unordered_2by2_lattice.txt', 2, include_analytical=True)
plotConvergence('../unordered_20by20_lattice.txt', 20)
plt.show()
