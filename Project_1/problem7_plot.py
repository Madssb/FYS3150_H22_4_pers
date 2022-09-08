import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('problem7_data.txt')
# print(data)
x = data[:,0]
u = data[:,1]

sol = 1 - (1 - np.exp(-10)) * x - np.exp(-10*x)


fig = plt.figure()
ax = fig.add_subplot()
ax.plot(x, u, 'r', label='Numerical', linestyle='--', linewidth=2)
ax.plot(x, sol, 'b', label='Analytical', linestyle=':', linewidth=2)
ax.set_xlabel('x', weight='bold', fontsize=16)

fig.tight_layout()
ax.legend(prop={'size': 16})
plt.show()
