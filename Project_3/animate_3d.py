import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

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

fig = plt.figure(figsize=(10, 8.5))
ax = plt.axes(projection='3d')

def animate(i):

    ax.clear()
    ax.grid(False)

    ax.plot3D(r1[:i*100+1, 0], r1[:i*100+1, 1], r1[:i*100+1, 2], color='red', lw=.25, ls='dashed')
    ax.scatter(r1[i*100, 0], r1[i*100, 1], r1[i*100, 2], color='red', lw=1, label='$P_1$')

    ax.plot3D(r2[:i*100+1, 0], r2[:i*100+1, 1], r2[:i*100+1, 2], color='royalblue', lw=.25, ls='dashed')
    ax.scatter(r2[i*100, 0], r2[i*100, 1], r2[i*100, 2], color='royalblue', lw=1, label='$P_2$')

    ax.set_title(f'Two particles in Penning trap\n$t={t[i*100]:.3f}\,\mu s$')
    ax.set_xlabel('x [$\mu s$]')
    ax.set_ylabel('y [$\mu s$]')
    ax.set_zlabel('z [$\mu s$]')
    ax.legend()

    ax.set_xlim3d([-40, 40])
    ax.set_ylim3d([-40, 40])
    ax.set_zlim3d([-30, 30])


ani = animation.FuncAnimation(fig, animate, N//100, interval=1)


# writergif = animation.PillowWriter(fps=N/100)
# ani.save('animation.gif', writer=writergif)

plt.show()
