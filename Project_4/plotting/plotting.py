# Code to visualize how estimates of <e> and <|m|> evolve with the number of
# MCMC cycles

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def analyticalValues():
    '''
    Return the analytical values of energy and magnetization per spin,
    specific heat cap. and susceptibility for 2x2 case
    '''

    # Expected values for <e>, <e^2>, <|m|> and <m^2>
    e = -2 * np.sinh(8) / (np.cosh(8) + 3)
    e_sqrd = 4 * np.cosh(8) / (np.cosh(8) + 3)
    m = 1. / 2 * (np.exp(8) + 2) / (np.cosh(8) +3)
    m_sqrd = 1. / 2 * (np.exp(8) + 1) / (np.cosh(8) + 3)
    C_V = 4 * (e_sqrd - e**2)
    X = 4 * (m_sqrd - m**2)

    return e, e_sqrd, m, m_sqrd, C_V, X

def plotConvergence(filename, L, T, cycles=1000000, include_analytical=None, parallel=None):
    '''
    Take a .txt file with computed values and plot the
    evolment with number of MCMC cycles
    '''

    nCycles = np.linspace(1, cycles, cycles)

    e, m, c, x, _, __ = np.loadtxt(filename, unpack=True)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 6.75), sharex=True)

    if parallel == True:

        fig.suptitle(f'{L}x{L} lattice at T={T:.1f} K\nfrom parallel coding')

    else:
        fig.suptitle(f'{L}x{L} lattice at T={T:.1f} K')

    ax1.set_title('Energy per spin')
    ax1.plot(nCycles, e, lw=1, color='red', label=r'$\overline{\epsilon}$')
    ax1.set_xscale('log')
    ax1.set_ylabel('[J]')

    ax2.set_title('Magnetization per spin')
    ax2.plot(nCycles, m, lw=1, color='royalblue', label='$\\overline{m}$')
    ax2.set_xscale('log')
    ax2.set_ylabel('M')

    ax3.set_title('Specific heat capacity')
    ax3.plot(nCycles, c, lw=1, color='black', label='$\\overline{C_V}$')
    ax3.set_ylabel('$C_V$')
    ax3.set_xscale('log')
    ax3.set_xlabel('No. of cycles')

    ax4.set_title('Susceptibility')
    ax4.plot(nCycles, x, lw=1, color='green', label='$\\overline{\chi}$')
    ax4.set_xscale('log')
    ax4.set_ylabel('$\chi$')
    ax4.set_xlabel('No. of cycles')

    if(include_analytical == True):

        e, e_sqrd, m, m_sqrd, C_V, X = analyticalValues()

        y = np.ones(len(nCycles))

        ax1.plot(nCycles, y * e, color='green', lw=1, ls='dashed', label='$<\\varepsilon>$')
        ax2.plot(nCycles, y * m, color='black', lw=1, ls='dashed', label='<|m|>')
        ax3.plot(nCycles, y * C_V, color='royalblue', lw=1, ls='dashed', label=r'$C_V^{analytical}$')
        ax4.plot(nCycles, y * X, color='red', lw=1, ls='dashed', label=r'$\chi^{analytical}$')

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()

    fig.tight_layout()

def plotOrderedUnordered(ordered_unordered, cycles=1000000, save=None):

    nCycles = np.linspace(1, cycles, cycles)

    fig, axes = plt.subplots(2, 2, figsize=(12, 6.75), sharex=True)
    fig.suptitle('20x20 lattice')

    T = [1, 2]
    colors = ['red', 'royalblue']
    labels = ['1', '2.4']
    ylabels = [r'Energy $[J]$', 'Magnetization']

    for i in range(2):

        filename = ordered_unordered + f'_L20_T{T[i]}.txt'

        for j in range(2):

            val = np.loadtxt(filename, usecols=j)

            axes[i, j].plot(nCycles, val, color=colors[i], lw=1, label=f'T={labels[i]}')
            axes[i, j].set_xlabel('Cycles')
            axes[i, j].set_ylabel(ylabels[j])
            axes[i, j].set_xscale('log')
            axes[i, j].legend()

    fig.subplots_adjust(hspace=0)


def histogram(ordered_unordered, L=20, burnIn=0, save=None):

    fig, axes = plt.subplots(2, 1, figsize=(12, 6.75))
    fig.suptitle('Probability density from ' + ordered_unordered + ' state')

    T = ['1', '2']

    for i in range(2):

        filename = ordered_unordered + '_L' + f'{L}' + '_T' + T[i] + '.txt'
        e = np.loadtxt(filename, usecols=4)
        eBurn = e[burnIn:]

        axes[i].hist(eBurn, bins='auto', histtype='stepfilled', density=True, color='royalblue', alpha=1)
        axes[i].set_title('T=' + T[i] + ' K')

    if save == True:

        plt.savefig('../figures/pdf_' + f'{L}.pdf')

    fig.tight_layout()


def burnIn(filename1, filename2, L, T1, T2, save=None):

    nCycles1, eps1, magn1, _, __ = np.loadtxt(filename1, unpack=True)
    nCycles2, eps2, magn2, _, __ = np.loadtxt(filename2, unpack=True)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, figsize=(8, 4.5))

    ax1.set_title(f'T={T1} K')
    ax1.plot(nCycles1, eps1, color='red', lw=1, label=r'$\varepsilon$')
    ax1.legend()

    ax2.set_title(f'T={T2} K')
    ax2.plot(nCycles2, eps2, color='red', lw=1, label=r'$\varepsilon$')
    ax2r = ax2.secondary_yaxis('right')
    ax2r.tick_params('y')
    ax2.set_yticks([])
    ax2.legend()

    ax3.plot(nCycles1, magn1, color='black', lw=1, label='<|m|>')
    ax3t = ax3.secondary_xaxis('top')
    ax3t.tick_params('x', direction='inout')
    ax3t.set_xticklabels([])
    ax3.legend()

    ax4.plot(nCycles2, magn2, color='black', lw=1, label='<|m|>')
    ax4t = ax4.secondary_xaxis('top')
    ax4t.tick_params('x', direction='inout')
    ax4t.set_xticklabels([])
    ax4r = ax4.secondary_yaxis('right')
    ax4r.tick_params('y')
    ax4.set_yticks([])
    ax4.legend()

    plt.subplots_adjust(hspace=0, wspace=0)


def plotParallel(filename, L, T, nruns, include_analytical=None, threads=None):

    fig, ax = plt.subplots(2, 1, figsize=(10.6, 6), sharex=True)
    title = f'{L}x{L} lattice at T={T:.1f} K'

    labels = np.array([r'$\epsilon$', 'm', r'$C_V$', r'$\chi$'])
    cycles = np.linspace(1, nruns, nruns)

    if threads != None:

        fig.suptitle(title + f'\nfrom {threads} threads')

        e1, m1, _, __, ___, ____ = np.loadtxt(filename, max_rows=nruns, unpack=True)
        e2, m2, _, __, ___, ____ = np.loadtxt(filename, skiprows=nruns, max_rows=nruns, unpack=True)
        e3, m3, _, __, ___, ____ = np.loadtxt(filename, skiprows=2 * nruns, max_rows=nruns, unpack=True)
        e4, m4, _, __, ___, ____ = np.loadtxt(filename, skiprows=3 * nruns, max_rows=nruns, unpack=True)

        e = np.array([e1, e2, e3, e4])
        m = np.array([m1, m2, m3, m4])

        ebase = e[0, -1]
        emin = e.min() - ebase
        emax = e.max() - ebase
        elim = [e.min() + abs(emin * .5), e.max() - abs(emax * .95)]

        mbase = m[0, -1]
        mmin = m.min() - mbase
        mmax = m.max() - mbase
        mlim = [m.min() + abs(mmin * .3), m.max() - abs(mmax * .3)]

        vals = np.array([e, m])

        for i in range(len(ax)):
            for j in range(len(e)):

                val = vals[i, j, :]

                ax[i].plot(cycles, val, lw=.75)

    else:

        fig.suptitle(title)

        e, m, _, __ = np.loadtxt(filename, unpack=True)
        elim = []
        mlim = []

        vals = np.array([e, m])

        for i in range(len(ax)):

            val = vals[i]

            ax[i].plot(cycles, val, lw=.75)

    ax[0].tick_params(axis='y', direction='in')

    ax[1].tick_params(axis='both', direction='in')
    ax1t = ax[1].secondary_xaxis('top')
    ax1t.tick_params(axis='x', direction='inout')
    ax1t.set_xticklabels([])

    # ax[0].set_ylim(elim)
    ax[0].set_ylabel('Energy per spin [J]')

    # ax[1].set_ylim(mlim)
    ax[1].set_xlabel('No. of cycles')
    ax[1].set_ylabel('Magnetization')

    plt.subplots_adjust(hspace=0, wspace=0)


def plotPhase(burn_in=None, save=None):

    N = 6
    L = [40, 60, 80, 100]
    T = np.linspace(2.1, 2.4, N)

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 6))
    ax = axes.flatten()

    lines = ['solid', 'dotted', 'dashed', 'dashdot']
    ylabels=[r'$<\epsilon>$', r'$<|m|>$', r'$C_V$', r'$\chi$']

    for i in range(4):

        filename = 'phase_L' + f'{L[i]}_v2.txt'

        for j in range(4):

            val = np.loadtxt(filename, usecols=j)

            ax[j].plot(T, val, color='black', lw=.75, ls=lines[i], label=f'{L[i]}x{L[i]}')
            ax[j].set_xlabel(r'Temperature $[T/k_B]$')
            ax[j].set_ylabel(ylabels[j])
            ax[j].legend()

    fig.tight_layout()

    if save == True:

        plt.savefig('../figures/phase_transition.pdf')


# plotConvergence('unordered_L2_T1.txt', 2, 1, include_analytical=True)
# plotConvergence('unordered_L20_T2.txt', 20, 2.4)
# plotOrderedUnordered('ordered')
# burnIn('../unordered_20by20_lattice_temp_1.txt', '../unordered_20by20_lattice_temp_2.txt', 20, 1, 2.4)
# plotParallel('../parallel_L40_v2.txt', 40, 1, 1000000, threads=4)
# plotParallel('../parallel_L20_T2.txt', 20, 2.5, 1000000, threads=4)
# plotParallel('../L20_T2.txt', 20, 2.5, 1000000)
# plotParallel('../test_1.txt', 60, 2.5, 100000, threads=4)
# plotPhase(save=True)


burnInIDX = 250000

# histogram('ordered', burnIn=burnInIDX)
# histogram('unordered', burnIn=burnInIDX)

plt.show()
