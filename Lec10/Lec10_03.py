"""Lec10_03: u_t + c(x,t) u_x = f(x,t) зі змінним коефіцієнтом (неявна центральна схема,
на кожному кроці СЛАР розв'язується немонотонною прогонкою)."""
import numpy as np
import matplotlib.pyplot as plt
from Lab08.bvp import m_progn


def c(x, t):
    return (t + 0.5) * (0.6 - x) ** 3


def f(x, t):
    return np.sin(x - t)


def u0(x):
    return x * (1 - x)


def g1(t):
    return t


def g2(t):
    return 2 * t


a, b = 0, 1
N = 101
x = np.linspace(a, b, N)
h = x[1] - x[0]
N1 = N - 1
tau, M = 0.01, 200
th = 0.5 * tau / h

A = np.zeros(N)
B = np.zeros(N)
C = -np.ones(N)
F = np.zeros(N)
t = 0.0
y = u0(x)
U0 = y.copy()
U = np.zeros((N, M + 1))
U[:, 0] = y
U1 = None
for j in range(1, M + 1):
    t += tau
    # коефіцієнти СЛАР
    F[0] = -g1(t)
    for k in range(1, N1):
        z = c(x[k], t) * th
        A[k] = -z
        B[k] = z
        F[k] = -(y[k] + tau * f(x[k], t))
    F[N - 1] = -g2(t)
    y, _ = m_progn(A, C, B, F)
    U[:, j] = y
    if 2 * j == M:
        U1 = y.copy()

# графіки розв'язку при t=0, M*tau/2, M*tau
t = M * tau
plt.plot(x, U0, 'b.-', x, U1, 'k.--', x, y, 'r*-')
plt.xlabel('x'); plt.ylabel('u')
plt.title("Наближений розв'язок")
plt.legend(['u(x,0)', f'u(x,{t / 2:g})', f'u(x,{t:g})'], loc='best')
plt.grid(True)
plt.show()

# поверхня u(x,t) з проекцією ліній рівня (аналог surfc)
tt = np.linspace(0, t, M + 1)
T, X = np.meshgrid(tt, x)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_surface(X, T, U, cmap='viridis', linewidth=0, antialiased=False)
ax.contour(X, T, U, zdir='z', offset=U.min(), cmap='viridis')
ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel('u')
ax.set_title('Наближений розв\'язок u(X,T)')
plt.show()
