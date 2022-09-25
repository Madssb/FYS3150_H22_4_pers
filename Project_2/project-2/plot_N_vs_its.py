# Plotting results from the function jacobi_eigensolver_multiple in souce.cpp
import numpy as np
import matplotlib.pyplot as plt

N, its = np.loadtxt("N_vs_iterations_DONT_MODIFY.txt", unpack=True)

plt.plot(N, its/N, color='royalblue')
plt.xticks(N)
plt.xlabel('N')
plt.ylabel('Iterations/N')
# plt.savefig('N_vs_iterations.pdf')
plt.show()

# Dividing iterations on N and printing to terminal
print(its/N)
"""
terminal> [ 1.          3.66666667  1.75        6.8         8.83333333 10.71428571
 12.         14.22222222 15.5        17.54545455 19.58333333 20.76923077
 22.85714286 24.8        26.5        28.29411765 29.72222222 30.42105263
 33.3       ]
"""
