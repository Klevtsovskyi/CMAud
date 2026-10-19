"""Pr10_1: модельна КЗ  U'' - p(x)U = -f(x),  p(x)=4+x,  x in [0,3].
f(x)=6+x(5+x)-(8+x)cos(x)^2, точний розв'язок u(x)=x+sin(x)^2."""
import numpy as np
import matplotlib.pyplot as plt
from bvp import m_grid2

a, b = 0, 3


def p(x):
    return 4 + x


def f(x):
    return 6 + x * (5 + x) - (8 + x) * np.cos(x) ** 2


# Варіанти крайових умов (лишаємо останній, як в оригіналі):
# D = [0, 1, 0, 0, 1, 3 + np.sin(b) ** 2]     # u(a)=-D3, u(b)=D6
# D = [1, 0, -1, 1, 0, 1 + np.sin(2 * 6)]     # u'(a)=-D3, u'(b)=D6
D = [2, -1, -2, 3, 2, 8.201583218752855]


def u(x):
    return x + np.sin(x) ** 2               # точний розв'язок


print('  h       ||y-u(X)||c')
k = 5
for n in (31, 301, 3001, 30001):
    X, Y = m_grid2(p, f, a, b, D, n)
    h = X[1] - X[0]
    U = u(X)
    m, M = U.min(), U.max()
    I = np.arange(0, n, k)
    k *= 10

    plt.figure()
    plt.plot(X, Y, 'b-', X[I], U[I], 'k*')
    plt.grid(True)
    plt.title(f'МРС O(h^2) h={h:6.4f}')
    de = 1e-1
    plt.axis([a - de, b + de, m - de, M + de])
    plt.xlabel('x'); plt.ylabel('y,u(x)')
    plt.legend(['y', 'u(x)'], loc='best')
    plt.show()

    E = np.abs(Y - U)
    plt.figure()
    plt.plot(X, E, 'r', linewidth=1.5)
    plt.grid(True)
    plt.title('Абсолютна похибка')
    plt.xlabel('x'); plt.ylabel('|y-u(x)|')
    print(f'{h:7.4f} {E.max():14.6e}')
    plt.show()
