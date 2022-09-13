import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot()
x0 = np.linspace(0, 1, 10000)
u = 1 - (1 - np.exp(-10))*x0 - np.exp(-10*x0)     # Analytic solution
axin = ax.inset_axes([0.185, 0.05, 0.25, 0.4])

n_list = [10**i for i in range(1, 7)]
for i, n in enumerate(n_list, start=1):
    data = np.loadtxt(f'problem7_data_n={n}.txt')   # Collect numerical data
    x = data[:, 0]
    v = data[:, 1]

    ax.plot(x, v, linestyle='--', label=f'Numeric, n = 10' +
            r'$^{' + f'{i}' + r'}$')
    axin.plot(x, v, linestyle='--')

ax.plot(x0, u, 'b', linestyle=':', label='Analytic solution')
axin.plot(x0, u, 'b', linestyle=':')
ax.set_xlabel('x', weight='bold', fontsize=24)
axin.set_xlim(0.243, 0.2433)
axin.set_ylim(0.6689, 0.6697)
axin.set_xticks([0.2430, 0.2431, 0.2432, 0.2433])
axin.set_yticks([0.6689, 0.6691, 0.6693, 0.6695, 0.6697])
ax.indicate_inset_zoom(axin, edgecolor="black", alpha=0.3)

ax.legend(prop={'size': 20})
fig.tight_layout()
plt.show()
