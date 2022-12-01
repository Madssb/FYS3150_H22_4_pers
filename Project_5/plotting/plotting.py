import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

plt.rc('font', size=14)
plt.rc('axes', labelsize=16)

def deviation(filename_in, T=2.5e-5, save=False):

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

    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(abs(S[0])))
    img = ax.imshow(abs(S[0,:,:]), extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("afmhot"), norm=norm)

    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label('Wave intensity', fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    time_txt = plt.text(0.95, 0.95, "t = {:.3f}".format(t_min), color="white",
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)


    def animation(i):
        # Normalize the colour scale to the current frame?
        # norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(abs(S[i,:,:])))
        # img.set_norm(norm)

        # Update z data
        img.set_data(abs(S[i,:,:]))

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3f}".format(current_time))

        return img

    anim = FuncAnimation(fig, animation, interval=10, frames=np.arange(0, len(S[0,:,:]), 2), repeat=True, blit=0)

    if save == True:

        writer = matplotlib.animation.PillowWriter(fps=30)
        anim.save('figures/' + filename_in + '.gif', writer=writer)

    else:

        plt.show()

deviation('problem_7', save=True)
deviation('problem_7_wo_potential', save=True)
# animate('problem_7', save=True)
# animate('problem_7')
# animate('problem_7_wo_potential', save=True)
