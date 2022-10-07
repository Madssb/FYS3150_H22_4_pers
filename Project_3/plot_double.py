import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

t, x1, x2, y1, y2 = np.loadtxt('data_test_double.txt', unpack=True)

fig, ax = plt.subplots()
line, = ax.plot(x1, y1)

def animate(i):

    line.set_data(x1[:i], y1[:i])

    return line,

ani = animation.FuncAnimation(fig, animate, len(x1), interval=200, blit=True)


# plt.title(f'Particle motion after $t={int(t[-1])}\,\mu s$')
#
# plt.plot(x1, y1, color='red', lw=1, label="$P_1$")
# plt.plot(x2, y2, color='royalblue', lw=1, label="$P_2$")
#
# plt.xlabel('x [$\mu m$]')
# plt.ylabel('y [$\mu m$]')
#
# plt.legend()
plt.show()
