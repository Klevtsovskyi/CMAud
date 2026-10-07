import numpy as np
import matplotlib.pyplot as plt
from inteq import fr2_mpn


def k(x, s):
    return x * s + 1               # K(x,s)


def f(x):
    return (x - 1 / 16) * x - 1 / 12


lamb, a, b = 0.25, 0, 1
n, eps, mi = 13, 1e-7, 5000
X, Y, n = fr2_mpn(lamb, a, b, 1, n, eps, mi, k, f)

plt.plot(X, Y, 'm', linewidth=2)
plt.title("Наближений розв'язок IP y")
plt.xlabel('x'); plt.ylabel('y'); plt.grid(True)
plt.show()

Z = X * X - Y                      # точний розв'язок u(x)=x^2
print(f'max|u-y| = {np.max(np.abs(Z)):.3e}')
plt.plot(X, Z, 'r', linewidth=2)
plt.title("Різниця розв'язків u-y")
plt.xlabel('x'); plt.ylabel('u-y'); plt.grid(True)
plt.show()
