import numpy as np
import matplotlib.pyplot as plt
from inteq import frur_mk


def k(x, t, u):
    y = np.exp(x - u)      # ядро K(x,t,u)
    return y, -y           # і dK/du


def f(x, u):
    return x - u, -1.0     # f(x,u) і df/du


a, b = 0, 0.35
n, eps, mi = 13, 1e-8, 500
Y0 = 0.05 * np.ones(n)     # початкове наближення
X, Y = frur_mk(a, b, 1, n, k, f, Y0, eps, mi, 2)

plt.plot(X, Y, 'k', linewidth=2)
plt.title("Наближений розв'язок IP")
plt.xlabel('x'); plt.ylabel('y'); plt.grid(True)
plt.show()
