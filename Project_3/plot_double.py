# File to plot figures
import numpy as np
import matplotlib.pyplot as plt

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

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 6.5))

fig.suptitle(f'Simulation for $t={np.ceil(t[-1])}\,\mu s$')

ax1.set_title('Motion')
ax1.plot(r1[:, 0], r1[:, 1], color='red', lw=.5, label="$P_1$")
ax1.plot(r2[:, 0], r2[:, 1], color='royalblue', lw=.5, label="$P_2$")

ax1.set_xlabel('x [$\mu m$]')
ax1.set_ylabel('y [$\mu m$]')
ax1.legend()

ax2.set_title('Phase space $(x,v_x)$')
ax2.plot(r1[:, 0], v1[:, 0], color='red', lw=.5, label="$P_1$")
ax2.plot(r2[:, 0], v2[:, 0], color='royalblue', lw=.5, label="$P_2$")

ax2.set_xlabel('x [$\mu m$]')
ax2.set_ylabel('$v_x$ [$\mu m/\mu s$]')
ax2.legend()

ax3.set_title('Phase space $(y,v_y)$')
ax3.plot(r1[:, 1], v1[:, 1], color='red', lw=.5, label="$P_1$")
ax3.plot(r2[:, 1], v2[:, 1], color='royalblue', lw=.5, label="$P_2$")

ax3.set_xlabel('y [$\mu m$]')
ax3.set_ylabel('$v_y$ [$\mu m/\mu s$]')
ax3.legend()

ax4.set_title('Phase space $(z,v_z)$')
ax4.plot(r1[:, 2], v1[:, 2], color='red', lw=.5, label="$P_1$")
ax4.plot(r2[:, 2], v2[:, 2], color='royalblue', lw=.5, label="$P_2$")

ax4.set_xlabel('z [$\mu m$]')
ax4.set_ylabel('$v_z$ [$\mu m/\mu s$]')
ax4.legend()



fig.tight_layout()
plt.show()
