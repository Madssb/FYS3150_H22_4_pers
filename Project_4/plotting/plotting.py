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
    m_sqrd = 1. / 2 * (np.exp(8) + 1) / (np.cosh(8) + 3)
    C_V = 4 * (e_sqrd - e**2)
    X = 4 * (m_sqrd - m**2)

    return e, e_sqrd, m, m_sqrd, C_V, X

def plotConvergence(filename, L, T, include_analytical=None, parallel=None):
    '''
    Take a .txt file with computed values and plot the
    evolment with number of MCMC cycles
    '''

    nCycles, eps, magn, C_V, X = np.loadtxt(filename, unpack=True)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10.6, 6), sharex=True)

    if parallel == True:

        fig.suptitle(f'{L}x{L} lattice at T={T:.1f} K\nfrom parallel coding')

    else:
        fig.suptitle(f'{L}x{L} lattice at T={T:.1f} K')

    ax1.set_title('Energy per spin')
    ax1.plot(nCycles, eps, lw=1, color='red', label=r'$\overline{\epsilon}$')
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

    plt.tight_layout()

def plotOrderedUnordered(orderedFilename, unorderedFilename, L, T, save=None):

    oCycles, oEps, oMagn, _, __ = np.loadtxt(orderedFilename, unpack=True)
    uCycles, uEps, uMagn, _, __ = np.loadtxt(unorderedFilename, unpack=True)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
    fig.suptitle(f'{L}x{L} lattice at T={T:.1f} K')

    ax1.plot(oCycles, oEps, color='red', lw=1, label='$\epsilon$ for ordered lattice')
    ax1.tick_params(axis='y', direction='in')
    ax1t = ax1.secondary_xaxis('top')
    ax1t.tick_params(axis='x', direction='in')
    ax1t.set_xticklabels([])
    ax1.legend()

    ax2.plot(uCycles, uEps, color='black', lw=1, label='$\epsilon$ for unordered lattice')
    ax2.tick_params(axis='both', direction='in')
    ax2t = ax2.secondary_xaxis('top')
    ax2t.tick_params(axis='x', direction='inout')
    ax2t.set_xticklabels([])
    ax2.legend()

    plt.subplots_adjust(hspace=0)


def histogram(filename, L, T, burnIn, save=None):

    _, eps, __, ___, ____ = np.loadtxt(filename, unpack=True)

    eps = eps[burnIn:]

    plt.figure()
    plt.hist(eps, bins='auto', density=True, color='royalblue', alpha=1, fill=True)


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

        e1, m1, _, __ = np.loadtxt(filename, max_rows=nruns, unpack=True)
        e2, m2, _, __ = np.loadtxt(filename, skiprows=nruns, max_rows=nruns, unpack=True)
        e3, m3, _, __ = np.loadtxt(filename, skiprows=2 * nruns, max_rows=nruns, unpack=True)
        e4, m4, _, __ = np.loadtxt(filename, skiprows=3 * nruns, max_rows=nruns, unpack=True)

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

    ax[0].set_ylim(elim)
    ax[0].set_ylabel('Energy per spin [J]')

    ax[1].set_ylim(mlim)
    ax[1].set_xlabel('No. of cycles')
    ax[1].set_ylabel('Magnetization')

    plt.subplots_adjust(hspace=0, wspace=0)

def plotPhase(filename):

    N = 6
    L = [40, 60, 80, 100]
    T = np.linspace(2.1, 2.4, N)

    c1, x1 = np.loadtxt(filename, max_rows=N, unpack=True)
    c2, x2 = np.loadtxt(filename, skiprows=N, max_rows=N, unpack=True)
    c3, x3 = np.loadtxt(filename, skiprows=2*N, max_rows=N, unpack=True)
    c4, x4 = np.loadtxt(filename, skiprows=3*N, max_rows=N, unpack=True)

    c = np.array([c1, c2, c3, c4])
    x = np.array([x1, x2, x3, x4])

    vals = np.array([c, x])

    fig, ax = plt.subplots(2, 1, figsize=(10, 5))

    for i in range(len(ax)):
        for j in range(len(c)):

            val = vals[i, j, :]

            ax[i].plot(T, val, lw=.75, label=f'L={L[j]}')
            ax[i].legend()

    fig.tight_layout()

# plotConvergence('../unordered_2by2_lattice_temp_1.txt', 2, 1,  include_analytical=True)
# plotOrderedUnordered('../unordered_20by20_lattice_temp_2.txt', '../ordered_20by20_lattice_temp_2.txt', 20, 2.4)
# burnIn('../unordered_20by20_lattice_temp_1.txt', '../unordered_20by20_lattice_temp_2.txt', 20, 1, 2.4)
# plotParallel('../parallel_L2_T1.txt', 2, 1, 1000000, threads=4)
plotParallel('../parallel_L20_T2.txt', 20, 2.5, 1000000, threads=4)
# plotParallel('../L20_T2.txt', 20, 2.5, 1000000)
# plotPhase('../test.txt')

T1BurnInIdx = 50000
T2BurnInIdx = 100000

# histogram('../unordered_20by20_lattice_temp_1.txt', 20, 1, T1BurnInIdx)
# histogram('../unordered_20by20_lattice_temp_2.txt', 20, 2.4, T2BurnInIdx)


plt.show()