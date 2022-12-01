import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def deviation(filename_in, T=2.5e-5, save=False):

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
    plt.xlabel('Time')
    plt.ylabel('Rel. error')
    plt.tight_layout()

    if save:

        plt.savefig(out_filename)

    plt.show()

def animate(filename_in, h=.005, dt=2.5e-5, save=None):

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

    filename = filename_in + '.bin'

    solution = pa.cx_cube()
    solution.load(filename)
    S = np.array(solution)

    t_list = np.array([0., .001, .002])

    fig, axes = plt.subplots(3, 3, figsize=(14, 7), sharex=True, sharey=True)
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(abs(S[0]))**2)
    norm_real = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(np.real(S[0])))
    norm_imag = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(np.imag(S[0])))

    for i in range(3):

        ax1 = axes[0, i]
        ax2 = axes[1, i]
        ax3 = axes[2, i]

        idx = int(t_list[i] / dt)
        S_n = S[idx, :, :]

        ax1.text(.95, .95, f't = {t_list[i]:.3f}', color='white', horizontalalignment='right', verticalalignment='top')
        ax2.text(.95, .95, f't = {t_list[i]:.3f}\n' + r'Re(u$_{ij}$)', color='white', horizontalalignment='right', verticalalignment='top')
        ax3.text(.95, .95, f't = {t_list[i]:.3f}\n' + r'Im(u$_{ij}$)', color='white', horizontalalignment='right', verticalalignment='top')
        img1 = ax1.imshow(abs(S_n)**2, extent=[0., 1., 0., 1.], cmap=plt.get_cmap("afmhot"), norm=norm)
        img2 = ax2.imshow(np.real(S_n), extent=[0., 1., 0., 1.], cmap=plt.get_cmap("afmhot"), norm=norm_real)
        img3 = ax3.imshow(np.imag(S_n), extent=[0., 1., 0., 1.], cmap=plt.get_cmap("afmhot"), norm=norm_imag)

        ax1.set_xlabel('x')
        ax2.set_xlabel('x')
        ax3.set_xlabel('x')

        ax1.set_ylabel('y')
        ax2.set_ylabel('y')
        ax3.set_ylabel('y')

        cbar1 =fig.colorbar(img1, ax=ax1)
        cbar2 =fig.colorbar(img1, ax=ax2)
        cbar3 =fig.colorbar(img1, ax=ax3)

    fig.tight_layout()
    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.75, 0.15, 0.02, 0.7])
    # fig.colorbar(img, cax=cbar_ax)
    # fig.subplots_adjust(hspace=0, wspace=0)

    if save:

        plt.savefig('figures/' + filename_in + '_real_imag.pdf')

    else:

        plt.show()

def probability(filename_in, out_filename_in, save=False):

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

    plt.figure(figsize=(8, 6))
    plt.plot(y, p, color='black', lw=1)
    plt.xlabel('y')
    plt.ylabel('Probability')

    if save:

        plt.savefig(out_filename)

    else:

        plt.show()

# deviation('problem_7', save=True)
# deviation('problem_7_wo_potential', save=True)
# animate('problem_7', save=True)
# animate('problem_7')
# animate('problem_7_wo_potential', save=True)
# animate('problem_8', save=True)
# colormaps_multiple('problem_8', save=True)
# colormaps_multiple('problem_8')
# probability('problem_8', 'prob_double_slit', save=True)
# probability('problem_9', 'prob_triple_slit', save=True)
probability('single_slit', 'prob_single_slit', save=True)
# animate('problem_9', save=True)
# animate('single_slit', save=True)
