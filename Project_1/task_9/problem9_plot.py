import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot()
x0 = np.linspace(0, 1, 10000)
u = 1 - (1 - np.exp(-10))*x0 - np.exp(-10*x0)     # Analytic solution
axin = ax.inset_axes([0.2, 0.1, 0.2, 0.2])

n_list = [10**i for i in range(1, 5)]
for n in n_list:
    data = np.loadtxt(f'problem9_data_n={n}.txt')   # Collect numerical data
    x = data[:,0]
    v = data[:,1]

    ax.plot(x, v, linestyle='--', label=f'Numeric, n = {n:.0e}')
    axin.plot(x, v, linestyle='--')

ax.plot(x0, u, 'b', linestyle=':', label='Analytic')
axin.plot(x0, u, 'b', linestyle=':')
ax.set_xlabel('x', weight='bold', fontsize=20)
axin.set_xlim(0.243, 0.2433)
axin.set_ylim(0.668, 0.677)
ax.indicate_inset_zoom(axin, edgecolor="black", alpha=0.3)

ax.legend(prop={'size': 14})
plt.show()
