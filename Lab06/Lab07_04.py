import numpy as np
import matplotlib.pyplot as plt
from inteq import vo2_mk

# а) модельна задача
def k_a(x, s):
    return s / (1 + x * x)


def f_a(x):
    return 1 / (1 + x * x)


lamb, a, b, n = -1, 0, 2, 51
X, Y, n = vo2_mk(lamb, a, b, n, k_a, f_a)
plt.plot(X, Y, 'm', linewidth=2)
plt.title("Наближений розв'язок IP y")
plt.xlabel('x'); plt.ylabel('y'); plt.grid(True)
plt.show()

Z = 1 / np.sqrt((1 + X * X) ** 3) - Y      # u - y
print(f'a) max|u-y| = {np.max(np.abs(Z)):.3e}')
plt.plot(X, Z, 'r')
plt.title("Різниця розв'язків u-y")
plt.xlabel('x'); plt.ylabel('u-y'); plt.grid(True)
plt.show()


# б)
def k_b(x, s):
    return x * s - 1


def f_b(x):
    return x


lamb, a, b = 1, 0, 2
X, Y, n = vo2_mk(lamb, a, b, n, k_b, f_b)
plt.plot(X, Y, 'k', linewidth=2)
plt.grid(True)
plt.title("Наближений розв'язок IP y")
plt.xlabel('x'); plt.ylabel('y')
plt.show()
