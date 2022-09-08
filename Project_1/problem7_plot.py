import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('problem7_data.txt')
x = data[:,0]
v = data[:,1]       # The approximated numerical solution

u = 1 - (1 - np.exp(-10)) * x - np.exp(-10*x)       # The exact solution


fig = plt.figure()
ax = fig.add_subplot()
ax.plot(x, v, 'r', label='Numerical', linestyle='--', linewidth=2)
ax.plot(x, u, 'b', label='Analytical', linestyle=':', linewidth=2)
ax.set_xlabel('x', weight='bold', fontsize=16)

fig.tight_layout()
ax.legend(prop={'size': 16})
plt.show()
