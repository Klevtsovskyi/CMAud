"""Pr10_2: модельна КЗ  U'' - 4U = -f(x),  f(x)=6+4x-8cos(x)^2,  x in [0,3].
Розв'язок u(x)=x+sin(x)^2. Порівняння схем O(h^2) і O(h^4)."""
import numpy as np
import matplotlib.pyplot as plt
from bvp import m_grid2, m_grid4

a, b = 0, 3
P = 4


def p(x):
    return 4


def f(x):
    return 6 + 4 * x - 8 * np.cos(x) ** 2


def u(x):
    return x + np.sin(x) ** 2               # точний розв'язок


def fd2f(x):
    """f(x), f''(x)"""
    return 6 + 4 * x - 8 * np.cos(x) ** 2, 16 * np.cos(2 * x)


def df1(x):
    """f'(x)"""
    return 4 + 8 * np.sin(2 * x)


# Варіанти крайових умов (I, II, III роду на кінцях):
# D = [0, 1, 0, 0, 1, u(b)]                 # I  р, I  р
# D = [0, 1, 0, 1, 0, 1 + np.sin(2*b)]      # I  р, II р
# D = [0, 1, 0, 3, 2, 8.201583218752855]    # I  р, III р
# D = [1, 0, -1, 0, 1, u(b)]                # II р, I  р
# D = [1, 0, -1, 1, 0, 1 + np.sin(2*b)]     # II р, II р
# D = [1, 0, -1, 3, 2, 8.201583218752855]   # II р, III р
# D = [2, 1, -2, 0, 1, u(b)]                # III р, I р
# D = [2, -1, -2, 1, 0, 1 + np.sin(2*b)]    # III р, II р
D = [2, -1, -2, 3, 2, 8.201583218752855]    # III р, III р

n, k = 301, 50
X, Y1 = m_grid2(p, f, a, b, D, n)
h = X[1] - X[0]
U = u(X)
m, M = U.min(), U.max()
I = np.arange(0, n, k)
de = 1e-1
lim = [a - de, b + de, m - de, M + de]

plt.figure()
plt.plot(X, Y1, 'b', linewidth=2)
plt.plot(X[I], U[I], 'k*')
plt.grid(True); plt.title('МСР O(h^2)'); plt.axis(lim)
plt.xlabel('x'); plt.ylabel('Y1,u(x)'); plt.legend(['Y1', 'u(x)'], loc='best')
plt.show()

E1 = np.abs(Y1 - U)
print('  h       ||Y-u(X)||c')
print(f'{h:7.4f} {E1.max():14.6e}')

X, Y2 = m_grid4(P, fd2f, df1, a, b, D, n)
plt.figure()
plt.plot(X, Y2, 'r', linewidth=2)
plt.plot(X[I], U[I], 'k*')
plt.grid(True); plt.title('МСР O(h^4)'); plt.axis(lim)
plt.xlabel('x'); plt.ylabel('Y2,u(x)'); plt.legend(['Y2', 'u(x)'], loc='best')
plt.show()

E2 = np.abs(Y2 - U)
print(f'{h:7.4f} {E2.max():14.6e}')

plt.figure()
plt.plot(X, E1, 'b', linewidth=1.5)
plt.plot(X, 1e5 * E2, 'r', linewidth=1.5)
plt.grid(True)
plt.title(r'Абсолютна похибка $|\psi|:=|Y-u(x)|$')
plt.xlabel('x'); plt.ylabel(r'$|\psi|$')
plt.legend([r'$|\psi|$ - O(h^2)', r'$10^5\cdot|\psi|$ - O(h^4)'], loc='best')
plt.show()

print('     U                         Y1                        Y2')
for ui, y1, y2 in zip(U[I], Y1[I], Y2[I]):
    print(f'{ui:.15e}  {y1:.15e}  {y2:.15e}')
