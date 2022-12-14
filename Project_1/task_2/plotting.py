"""
Displays data generated by "problem2.cpp".
"""
import numpy as np
import matplotlib.pyplot as plt

sol_vec = np.zeros((101, 2))  # Array to contain both x and u-values

with open('analytic_sol.data', 'r') as infile:

    for i in range(3):  # Reads the first 3 lines that are not important
        infile.readline()

    i = 0   # Counter
    for line in infile:
        data = line.rstrip().split('\t')
        sol_vec[i, 0], sol_vec[i, 1] = [float(num) for num in data]
        i += 1

# Plotting
cm = 1/2.54
fig, ax = plt.subplots(figsize=(10*cm, 10*cm))
ax.plot(sol_vec[:, 0], sol_vec[:, 1], label='$u(x)$')
ax.set_title('Analytic solution')
ax.set_xlabel('$x$')
ax.set_ylabel('$u(x)$')
ax.legend()
plt.savefig('problem_2.pdf')
#plt.show()
