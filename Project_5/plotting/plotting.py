import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def animate(filename, h=.005, dt=2.5e-5):

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
    img = ax.imshow(abs(S[0,:,:]), extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("z(x,y,t)", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white",
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)


    def animation(i):
        # Normalize the colour scale to the current frame?
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(abs(S[i,:,:])))
        img.set_norm(norm)

        # Update z data
        img.set_data(abs(S[i,:,:]))

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(S[0,:,:]), 2), repeat=True, blit=0)

    plt.show()


# solution = pa.cx_cube()
# solution.load('problem_7.bin')
# S = np.array(solution)
# # S = S.reshape((S.shape[1], S.shape[2], S.shape[0]))
# print(S.shape)
# plt.imshow(np.abs(S[0]))
animate('problem_7.bin')
