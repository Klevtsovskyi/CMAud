import numpy as np
import matplotlib.pyplot as plt
from inteq import fr2_mcol, fr2_mk


def k(x, s):
    return 1 / (x + 1)                      # K(x,s)


def f(x):
    return x * x * (1 - x) - 1 / (x + 1)    # f(x)


lamb, a, b, n = 12, 0, 1, 101
X = np.linspace(a, b, n)
x, YY = fr2_mcol(lamb, X, k, f)
X, Y, n = fr2_mk(lamb, a, b, 1, n, k, f)

plt.plot(X, Y, 'm-', x, YY, 'b:')
plt.title("Наближений розв'язок IP y")
plt.legend(['м.КФ', 'м.колок'], loc='best')
plt.xlabel('x'); plt.ylabel('y'); plt.grid(True)
plt.show()


def u(x):
    return x * x * (1 - x)                  # точний розв'язок


print(f'max|u-y|: м.КФ = {np.max(np.abs(u(X) - Y)):.3e}, '
      f'м.колок = {np.max(np.abs(u(x) - YY)):.3e}')

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(X, u(X) - Y, 'r', linewidth=2)
plt.title("м.КФ: різниця розв'язків u-y")
plt.xlabel('x'); plt.ylabel('u-y'); plt.grid(True)
plt.subplot(2, 1, 2)
plt.plot(x, u(x) - YY, 'r', linewidth=2)
plt.title("м.Колокацій: різниця розв'язків u-y")
plt.xlabel('x'); plt.ylabel('u-y'); plt.grid(True)
plt.tight_layout()
plt.show()
