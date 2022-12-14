'''
This program contains plotting functions to visualise the wave packet
simulations.
'''

import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.signal import find_peaks


def deviation(filename_in, T=0.008, save=False):
    '''
    Computing the total probability relative to one at each time step
    '''

    plt.rc('font', size=14)
    plt.rc('axes', labelsize=16)

    filename = filename_in + '.bin'
    out_filename = 'figures/' + filename_in + '_dev.pdf'

    solution = pa.cx_cube()
    solution.load(filename)
    S = np.array(solution)

    N = len(S)
    t = np.linspace(0, T, N)
    devs = np.zeros(N)

    for n in range(N):

        S_n = S[n, :, :]
        dev = np.real(np.sum(np.conj(S_n) * S_n))
        devs[n] = dev

    abs_errs = abs(1 - devs)

    fig = plt.figure(figsize=(8, 4))
    plt.plot(t, abs_errs, color='black', lw=1)
    plt.xlabel('Time [s]')
    plt.ylabel('Rel. error')
    plt.tight_layout()

    print(f'Init. err.: {abs_errs[0]:.3e}\nFinal err.:{abs_errs[-1]:.3e}')

    if save:

        plt.savefig(out_filename)

    plt.show()


def animate(filename_in, h=.005, dt=2.5e-5, save=None):
    '''
    Produce animated .gif files of the wave packet traversing
    a potential barrier
    '''

    plt.rc('font', size=14)
    plt.rc('axes', labelsize=16)

    filename = filename_in + '.bin'

    solution = pa.cx_cube()
    solution.load(filename)
    S = np.array(solution)

    X = np.arange(0, 1 + h, h)
    Y = np.arange(0, 1 + h, h)
    x, y = np.meshgrid(X, Y, sparse=True)

    t = np.arange(0, 1 + dt, dt)

    fontsize = 12
    t_min = t[0]
    x_min, x_max = X[0], X[-1]
    y_min, y_max = Y[0], Y[-1]

    fig = plt.figure()
    ax = plt.gca()

    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(abs(S[0]))**2)
    img = ax.imshow(abs(S[0,:,:])**2, extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("afmhot"), norm=norm)

    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label('Probability', fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    time_txt = plt.text(0.95, 0.95, "t = {:.3f}".format(t_min), color="white",
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)


    def animation(i):

        img.set_data(abs(S[i,:,:])**2)

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3f}".format(current_time))

        return img

    anim = FuncAnimation(fig, animation, interval=10, frames=np.arange(0, len(S[:,0,0]), 2), repeat=True, blit=0)

    if save == True:

        writer = matplotlib.animation.PillowWriter(fps=30)
        anim.save('figures/' + filename_in + '.gif', writer=writer)

    else:

        plt.show()


def colormaps_multiple(filename_in, dt=2.5e-5, save=False):
    '''
    Will plot 9 subplots of the wave packet: 3 for the full probability,
    3 for the real part of u, and 3 for the imaginary part.
    '''
    filename = filename_in + '.bin'

    solution = pa.cx_cube()
    solution.load(filename)
    S = np.array(solution)

    S_cx = abs(S)**2
    S_real = np.real(S)
    S_imag = np.imag(S)

    S_list = np.array([S_cx, S_real, S_imag])
    x_ticks = np.linspace(.2, .8, 4)
    x_labels = ['0.2', '0.4', '0.6', '0.8']
    S_names = [r'$p_{ij}^n$', r'Re{$u_{ij}^n$}', r'Im{$u_{ij}^n$}']

    t_list = np.array([0., .001, .002])

    plt.rc('xtick', top=True)
    plt.rc('xtick', labeltop=True)
    plt.rc('xtick', bottom=False)
    plt.rc('xtick', labelbottom=False)
    fig, axes = plt.subplots(3, 3, figsize=(12, 9), sharex=True, sharey=True)

    for i in range(3):

        ax = axes[i, :]
        S_n = S_list[i]
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_n[0]))

        for j in range(3):

            idx = int(t_list[j] / dt)
            S_n_idx = S_n[idx, :, :]

            ax[j].text(.95, .95, f't = {t_list[j]:.3f}\n{S_names[i]}', color='white', horizontalalignment='right', verticalalignment='top')
            im = ax[j].imshow(S_n_idx, extent=[0., 1., 0., 1.], aspect='auto', cmap=plt.get_cmap("afmhot"), norm=norm)
            ax[j].set_xticks(x_ticks)
            ax[j].set_xticklabels(x_labels)

        cax = fig.add_axes([0.125, .65 - (.2765 * i), .775, .01])
        fig.colorbar(im, cax=cax, orientation='horizontal')


    fig.subplots_adjust(wspace=.007, hspace=.275)
    fig.text(.5, .04, 'Probability', ha='center')
    fig.text(.07, .5, 'Position in y direction', va='center', rotation='vertical')
    fig.text(.5, .93, 'Position in x direction', ha='center')


    if save:

        plt.savefig('figures/' + filename_in + '_real_imag.pdf')

    else:

        plt.show()


def colormaps(filename_in, dt=2.5e-5, save=False):
    '''
    Will produce 3 subplots, so have to be run 3 times to get full
    probability, real part of u and imaginary part of u
    '''

    filename = filename_in + '.bin'

    solution = pa.cx_cube()
    solution.load(filename)
    S = np.array(solution)

    S_cx = abs(S)**2
    S_real = np.real(S)
    S_imag = np.imag(S)

    S_list = np.array([S_cx, S_real, S_imag])
    t_list = np.array([0., .001, .002])
    savefig_names = ['complex', 'real', 'imag']
    x_ticks = np.linspace(.2, .8, 4)
    x_labels = ['0.2', '0.4', '0.6', '0.8']

    for n in range(len(S_list)):

        plt.rc('xtick', top=True)
        plt.rc('xtick', labeltop=True)
        plt.rc('xtick', bottom=False)
        plt.rc('xtick', labelbottom=False)

        S_n = S_list[n]

        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(S_n[0]))

        fig, axes = plt.subplots(1, 3, figsize=(12, 6), sharey=True)

        for i in range(len(axes)):

            ax = axes[i]

            idx = int(t_list[i] / dt)
            S_n_idx = S_n[idx, :, :]

            ax.text(.95, .95, f't = {t_list[i]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
            im = ax.imshow(S_n_idx, extent=[0., 1., 0., 1.], cmap=plt.get_cmap("afmhot"), norm=norm)
            ax.set_xticks(x_ticks)
            ax.set_xticklabels(x_labels)
            ax.set_xlabel('x', labelpad=-250)

        axes[0].set_ylabel('y')

        fig.subplots_adjust(wspace=.005)
        c_ax = fig.add_axes([0.125, .2, .775, .033])
        fig.colorbar(im, cax=c_ax, orientation='horizontal', label='Probability')

        if save:

            plt.savefig('figures/' + savefig_names[n] + '.pdf')

        else:

            continue

    if save:

        pass

    else:

        plt.show()


def probability(filename_in, out_filename_in, save=False):
    '''
    Computing and plotting the 1D probability function at time t=0.002 s,
    measured at x=0.8
    '''

    plt.rc('font', size=14)
    plt.rc('axes', labelsize=16)

    filename = filename_in + '.bin'
    out_filename = 'figures/' + out_filename_in + '.pdf'

    solution = pa.cx_cube()
    solution.load(filename)
    S = np.array(solution)

    t_idx = int(.002 / 2.5e-5)
    x_idx = int(.8 / .005)

    P = np.conj(S) * S
    p = P[t_idx, :, x_idx]
    p /= np.sum(p)
    y = np.linspace(0, 1, len(S[0]))

    peaks, _ = find_peaks(p, height=.005)

    plt.figure(figsize=(8, 6))
    plt.plot(y, p, color='black', lw=1)

    plt.xlabel('y')
    plt.ylabel('Probability')

    if save:

        plt.savefig(out_filename)

    else:

        print(f'Peaks for ' + out_filename_in + ':')

        for i in range(len(peaks)):
            print(f'p({y[peaks[i]]:.3f}) = {np.real(p[peaks[i]]):.3f}')

        plt.show()


deviation('double_slit_long', save=True)
deviation('no_potential', save=True)
animate('double_slit_long', save=True)
animate('no_potential', save=True)
animate('double_slit_short', save=True)
colormaps_multiple('double_slit_short', save=True)
probability('double_slit_short', 'prob_double_slit', save=True)
probability('triple_slit', 'prob_triple_slit', save=True)
probability('single_slit', 'prob_single_slit', save=True)
animate('triple_slit', save=True)
animate('single_slit', save=True)
colormaps('double_slit_short', save=True)
