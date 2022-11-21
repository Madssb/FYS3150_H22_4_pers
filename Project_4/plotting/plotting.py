'''
All of the programs in this file will need to have the .txt files in the same
folder as the program itself. They will use them to produce different plots
depending on the function. All plots will be saved to a folder called figures/
in the same folder as this program.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats 

def analyticalValues(T):
    '''
    Return the analytical values of energy and magnetization per spin,
    specific heat cap. and susceptibility for 2x2 case
    '''

    # Expected values for <e>, <e^2>, <|m|> and <m^2>
    e = -2 * np.sinh(8 / T) / (np.cosh(8 / T) + 3)
    e_sqrd = 4 * np.cosh(8 / T) / (np.cosh(8 / T) + 3)
    m = 1. / 2 * (np.exp(8 / T) + 2) / (np.cosh(8 / T) +3)
    m_sqrd = 1. / 2 * (np.exp(8 / T) + 1) / (np.cosh(8 / T) + 3)
    C_V = 4 * (e_sqrd - e**2)
    X = 4 * (m_sqrd - m**2)

    return e, e_sqrd, m, m_sqrd, C_V, X


def plotConvergence(filename, L, T, cycles=1000000, include_analytical=None, parallel=None, save=None):
    '''
    Take a .txt file with computed values and plot the
    evolment with number of MCMC cycles
    '''

    nCycles = np.linspace(1, cycles, cycles)

    e, m, c, x, _, __ = np.loadtxt(filename, unpack=True)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 4), sharex=True)

    ax1.plot(nCycles, e, lw=1, color='red', label=r'$\overline{\epsilon}$')
    ax1.set_xscale('log')
    ax1.set_ylabel('[J]')

    ax2.plot(nCycles, m, lw=1, color='royalblue', label='$\\overline{m}$')
    ax2.set_xscale('log')
    ax2.tick_params(axis='y', which='both', labelleft=False, labelright=True, right=True, left=False)
    ax2.set_ylabel(r'$<|m|>$', labelpad=-40 * 8 - 8)

    ax3.plot(nCycles, c, lw=1, color='black', label='$\\overline{C_V}$')
    ax3.set_ylabel(r'$C_V$')
    ax3.set_xscale('log')
    ax3.set_xlabel('No. of cycles')

    ax4.plot(nCycles, x, lw=1, color='green', label='$\\overline{\chi}$')
    ax4.set_xscale('log')
    ax4.tick_params(axis='y', which='both', labelleft=False, labelright=True, right=True, left=False)
    ax4.set_ylabel(r'$\chi$', labelpad= -40 * 8 - 8)
    ax4.set_xlabel('No. of cycles')

    if(include_analytical == True):

        e, e_sqrd, m, m_sqrd, C_V, X = analyticalValues(T)

        y = np.ones(len(nCycles))

        ax1.plot(nCycles, y * e, color='green', lw=1, ls='dashed', label='$<\\varepsilon>$')
        ax2.plot(nCycles, y * m, color='black', lw=1, ls='dashed', label='<|m|>')
        ax3.plot(nCycles, y * C_V, color='royalblue', lw=1, ls='dashed', label=r'$C_V^{analytical}$')
        ax4.plot(nCycles, y * X, color='red', lw=1, ls='dashed', label=r'$\chi^{analytical}$')

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()

    plt.subplots_adjust(hspace=0, wspace=0)

    if save == True:

        plt.savefig(f'../figures/convergence_L{L}_T{int(T)}.pdf')


def plotOrderedUnordered(ordered_unordered, cycles=1000000, save=None):
    '''
    This function plots the observables computed on two systems; one where all
    of the spins has been set to +1 and one where they are randomly initialized
    '''

    nCycles = np.linspace(1, cycles, cycles)

    fig, axes = plt.subplots(2, 2, figsize=(10, 4), sharex=True)

    T = [1, 2]
    colors = ['red', 'royalblue']
    labels = ['1', '2.4']
    ylabels = [r'$[J]$', r'$<|m|>$']

    for i in range(2):

        filename = ordered_unordered + f'_L20_T{T[i]}.txt'

        for j in range(2):

            val = np.loadtxt(filename, usecols=j)

            axes[i, j].plot(nCycles, val, color=colors[i], lw=1, label=f'T={T[i]}\n{ordered_unordered}')
            axes[i, j].set_xlabel('No. of cycles')

            if j == 1 or j == 3:

                axes[i, j].tick_params(axis='y', which='both', labelleft=False, labelright=True, right=True, left=False)
                axes[i, j].set_ylabel(ylabels[j], labelpad= -40 * 8 - 13)

            else:

                axes[i, j].set_ylabel(ylabels[j])
                axes[i, j].set_xscale('log')

            axes[i, j].legend()

    fig.subplots_adjust(hspace=0, wspace=0)

    if save == True:

        plt.savefig('../figures/' + ordered_unordered + '_L20.pdf')


def histogram(ordered_unordered, L=20, burnIn=0, save=None):
    '''
    Creating normalized histograms from either ordered or unordered systems
    '''

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 6))

    T = ['1', '2']

    for i in range(2):

        filename = ordered_unordered + '_L' + f'{L}' + '_T' + T[i] + '.txt'
        e = np.loadtxt(filename, usecols=4)
        eBurn = e[burnIn:]

        axes[i].hist(eBurn, bins='auto', histtype='stepfilled', density=True, color='royalblue', alpha=1)
        axes[i].set_xlabel(r'Temperature [J/k$_B$]', fontsize=18)

    axes[0].set_ylabel('Frequency', fontsize=18)
    axes[0].tick_params(axis='both', which='both', labelsize=18)

    axes[1].tick_params(axis='y', which='both', labelleft=False, labelright=True, right=True, left=False, labelsize=18)
    axes[1].tick_params(axis='x', which='both', labelsize=18)

    plt.subplots_adjust(wspace=0)

    if save == True:

        plt.savefig('../figures/pdf_' + f'{L}.pdf')


def plotPhase(wide_narrow, start_temp=2.1, stop_temp=2.4, n_points=6, burn_in=None, save=None, include_160=None):
    '''
    This function plots <e>, <|m|>, C_V and X in hopes of seeing signs of
    phase transition
    '''

    if include_160 == True:

        L = [60, 80, 100, 160]

    else:

        L = [40, 60, 80, 100]

    T = np.linspace(start_temp, stop_temp, n_points)

    fig, axes = plt.subplots(2, 2, figsize=(10, 4.5), sharex=True)
    ax = axes.flatten()

    lines = ['solid', 'dotted', 'dashed', 'dashdot']
    ylabels=[r'$\langle\epsilon\rangle$', r'$\langle|m|\rangle$', r'$C_V$', r'$\chi$']

    for i in range(4):

        filename = wide_narrow + '_phase_L' + f'{L[i]}.txt'

        for j in range(4):

            val = np.loadtxt(filename, usecols=j)

            ax[j].plot(T, val, color='black', lw=.75, ls=lines[i], label=f'{L[i]}x{L[i]}')
            ax[j].set_xlabel(r'Temperature [J/k$_B$]')

            if i == 1 or i == 3:

                ax[i].tick_params(axis='y', which='both', labelleft=False, labelright=True, right=True, left=False)
                ax[i].set_ylabel(ylabels[j], labelpad= -40 * 8 - 8)

            else:

                ax[j].set_ylabel(ylabels[j])

    ax[0].legend()

    plt.subplots_adjust(hspace=0, wspace=0)

    if save == True:

        plt.savefig(f'../figures/{wide_narrow}_phase_transition.pdf')


def criticalTemperature(wide_narrow, save=None):
    '''
    Estimate the critical temperature of a system of infinate
    lattice size.
    '''

    if wide_narrow == 'wide':

        temps = np.linspace(2.1, 2.4, 10)

    else:

        temps = np.linspace(2.2, 2.35, 10)

    fig, ax = plt.subplots(figsize=(7, 7))

    L = np.array([40, 60, 80, 100])
    tc = np.zeros(len(L))

    for i in range(len(L)):

        filename = wide_narrow + f'_phase_L{L[i]}.txt'

        c = np.loadtxt(filename, usecols=2)

        idx = np.where(c == np.max(c))[0][0]
        tc[i] = temps[idx]

    linreg = reg(1 / L, tc)
    a = linreg.slope
    b = linreg.intercept

    ax.plot(L, a / L + b, color='black', lw=1)




# print(analyticalValues(1)[4:])
# print(analyticalValues(2.4)[4:])

# plotConvergence('unordered_L2_T1.txt', 2, 1, include_analytical=True, save=True)
# plotConvergence('unordered_L20_T2.txt', 20, 2.4, save=True)

# plotOrderedUnordered('ordered', save=True)
# plotOrderedUnordered('unordered', save=True)

burnInIDX = 250000

# histogram('ordered', burnIn=burnInIDX)
# histogram('unordered', burnIn=burnInIDX, save=True)

# plotParallel('../parallel_L40_v2.txt', 40, 1, 1000000, threads=4)
# plotParallel('../parallel_L20_T2.txt', 20, 2.5, 1000000, threads=4)
# plotParallel('../L20_T2.txt', 20, 2.5, 1000000)
# plotParallel('../test_1.txt', 60, 2.5, 100000, threads=4)

# plotPhase('wide', n_points=10)
# plotPhase('narrow', n_points=10, start_temp=2.2, stop_temp=2.35)
# plotPhase('wide', n_points=10, save=True)
# plotPhase('narrow', n_points=10, start_temp=2.2, stop_temp=2.35, save=True)

criticalTemperature('wide')

plt.show()
