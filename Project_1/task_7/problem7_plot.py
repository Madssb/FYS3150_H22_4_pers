import numpy as np
import matplotlib.pyplot as plt

data_n10 = np.loadtxt('problem7_data_n=10.txt')
data_n100 = np.loadtxt('problem7_data_n=100.txt')
data_n1000 = np.loadtxt('problem7_data_n=1000.txt')

x_n10 = data_n10[:,0]
v_n10 = data_n10[:,1]       # The approximated numerical solution

x_n100 = data_n100[:,0]
v_n100 = data_n100[:,1]       # The approximated numerical solution

x_n1000 = data_n1000[:,0]
v_n1000 = data_n1000[:,1]       # The approximated numerical solution

u_n10 = 1 - (1 - np.exp(-10)) * x_n10 - np.exp(-10*x_n10)       # The exact solution
u_n100 = 1 - (1 - np.exp(-10)) * x_n100 - np.exp(-10*x_n100)       # The exact solution
u_n1000 = 1 - (1 - np.exp(-10)) * x_n1000 - np.exp(-10*x_n1000)       # The exact solution


fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

ax1.plot(x_n10, v_n10, 'r', label='Numerical', linestyle='--', linewidth=2)
ax1.plot(x_n10, u_n10, 'b', label='Analytical', linestyle=':', linewidth=2)
ax1.set_xlabel('x', weight='bold', fontsize=16)
ax1.set_title('n=10', weight='bold', fontsize=20)

ax2.plot(x_n100, v_n100, 'r', label='Numerical', linestyle='--', linewidth=2)
ax2.plot(x_n100, u_n100, 'b', label='Analytical', linestyle=':', linewidth=2)
ax2.set_xlabel('x', weight='bold', fontsize=16)
ax2.set_title('n=100', weight='bold', fontsize=20)

ax3.plot(x_n1000, v_n1000, 'r', label='Numerical', linestyle='--', linewidth=2)
ax3.plot(x_n1000, u_n1000, 'b', label='Analytical', linestyle=':', linewidth=2)
ax3.set_xlabel('x', weight='bold', fontsize=16)
ax3.set_title('n=1000', weight='bold', fontsize=20)

fig.tight_layout()
ax1.legend(prop={'size': 16})
ax2.legend(prop={'size': 16})
ax3.legend(prop={'size': 16})
plt.show()
