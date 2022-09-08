import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

n = 100
h = 1/(n)

a, b, c = (-1., 2., -1.)

x = np.linspace(0, 1, n+1)

u = 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)

def u_sol(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)

def g(x):
    return h**2 * (u_sol(x) - a/b * u_sol(x))

b_tilde = np.zeros(n+1)
g_tilde = np.zeros_like(b_tilde)
v = np.zeros_like(b_tilde)

b_tilde[0] = b
g_tilde[0] = h**2 * u[0]

for i in range(1, n+1):
    g_i = g(x[i])
    b_tilde[i] = b - a/b_tilde[i-1] * c
    g_tilde[i] = g_i - a/b_tilde[i-1] * g_tilde[i-1]

v[n] = g_tilde[n] / b_tilde[n]

for i in range(n, 0, -1):
    v[i-1] = (g_tilde[i-1] - c * v[i]) / b_tilde[i-1]

plt.plot(x, u, label='$u(x)$')
plt.plot(x, v, label='$v(x)$')

plt.legend()
plt.show()
