import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

_t, x1, x2, y1, y2, z1, z2, vx1, vx2, vy1, vy2, vz1, vz2 = np.loadtxt('data_test_double.txt', unpack=True)

t = _t[::2]
N = len(t)

r1 = np.zeros((N, 3))
r1[:, 0] = x1[::2]
r1[:, 1] = y1[::2]
r1[:, 2] = z1[::2]

v1 = np.zeros((N, 3))
v1[:, 0] = vx1[::2]
v1[:, 1] = vy1[::2]
v1[:, 2] = vz1[::2]

r2 = np.zeros((N, 3))
r2[:, 0] = x2[::2]
r2[:, 1] = y2[::2]
r2[:, 2] = z2[::2]

v2 = np.zeros((N, 3))
v2[:, 0] = vx2[::2]
v2[:, 1] = vy2[::2]
v2[:, 2] = vz2[::2]

x_min = np.min(r1[:, 0]) - 1
x_max = np.max(r1[:, 0]) + 1
y_min = np.min(r1[:, 1]) - 1
y_max = np.max(r1[:, 1]) + 1

fig, ax = plt.subplots()

def animate(i):

    ax.clear()

    ax.plot(r1[:i*20+1, 0], r1[:i*20+1, 1], ls='dashed', color='red', lw=.5)
    ax.plot(r1[i*20, 0], r1[i*20, 1], 'o', color='red', label='$P_1$')    # ms: markersize

    ax.plot(r2[:i*20+1, 0], r2[:i*20+1, 1], ls='dashed', color='royalblue', lw=.5)
    ax.plot(r2[i*20, 0], r2[i*20, 1], 'o', color='royalblue', label='$P_2$')

    ax.set_title(f'Two particles in Penning trap\n$t={t[i]:.3f}\,\mu s$')
    ax.set_xlabel('x [$\mu m$]')
    ax.set_ylabel('y [$\mu m$]')

    ax.set_xlim([-40, 40])
    ax.set_ylim([-40, 40])
    ax.legend()

ani = animation.FuncAnimation(fig, animate, N//20, interval=1)

plt.show()
