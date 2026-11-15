"""Le14_1: задача Діріхле для рівняння Пуассона в квадраті. Порівняння ітераційних методів
Якобі, Зейделя і верхньої релаксації (кількість ітерацій залежно від кроку сітки)."""
import time
import numpy as np
import matplotlib.pyplot as plt
from elliptic import im_jzr, view_u

a = b = 1
PAUSE = 2          # тривалість показу кожного графіка, с (0 - чекати закриття вікна)


def u(x, y):                       # точний розв'язок
    return 1 + 15 * (x * (a - x) * y * (b - y))


def f(x, y):                       # права частина
    return 30 * (x * (a - x) + y * (b - y))


def muh(x):
    return 1.0, 1.0


def muv(y):
    return 1.0, 1.0


ep, itmax = 1e-8, 30000
Lep = np.log10(1 / ep)

N = np.array([5, 11, 21, 41, 81])
Kit = np.zeros((len(N), 3))


def row(mc, it, n0, t, err, ro=None):
    s = f'{mc:>12} {it:5d} {n0:5d}{t:9.3f} {err:13.5e}'
    print(s if ro is None else f'{s} {ro:7.5f}')


for n, Nx in enumerate(N):
    Ny = Nx
    X = np.linspace(0, a, Nx)
    Y = np.linspace(0, b, Ny)
    Hx = X[1] - X[0]
    n0e = 2 * Lep / (np.pi * Hx)
    UU = u(X[None, :], Y[:, None])

    print(f'Nx=Ny={Nx:4d} Hx=Hy={Hx:7.5f}')
    print('Ітерац.метод Ітер. ІтерT Час(сек) max(max|U-u|) ro')

    t0 = time.perf_counter()
    U, it = im_jzr(0, f, muh, muv, X, Y, 0, ep, itmax)
    t = time.perf_counter() - t0
    Kit[n, 0] = it
    row('м.Якобі', it, int(n0e / (np.pi * Hx)), t, np.max(np.abs(U - UU)))
    view_u(PAUSE, X, Y, U, 15, UU)

    t0 = time.perf_counter()
    U, it = im_jzr(1, f, muh, muv, X, Y, 1, ep, itmax)
    t = time.perf_counter() - t0
    Kit[n, 1] = it
    row('м.Зейделя', it, int(n0e / (2 * np.pi * Hx)), t, np.max(np.abs(U - UU)))
    view_u(PAUSE, X, Y, U, 15, UU)

    ro = 1 / (0.5 + np.sin(0.5 * np.pi * Hx))        # оптимальний параметр
    t0 = time.perf_counter()
    U, it = im_jzr(2, f, muh, muv, X, Y, ro, ep, itmax)
    t = time.perf_counter() - t0
    Kit[n, 2] = it
    row('м.релаксації', it, int(n0e), t, np.max(np.abs(U - UU)), ro)
    view_u(PAUSE, X, Y, U, 15, UU)

# залежність кількості ітерацій від h^{-1}
plt.close('all')
Nm = N - 1
M = np.linspace(Nm[0], Nm[-1], 101)
r = 2 * Lep / np.pi * M
J = r / np.pi * M
Z = 0.5 * J
plt.plot(Nm, Kit[:, 0], 'ks', linewidth=2)
plt.plot(Nm, Kit[:, 1], 'ko', linewidth=2)
plt.plot(Nm, Kit[:, 2], 'b*', linewidth=2)
plt.plot(M, J, 'r-.', linewidth=2)
plt.plot(M, Z, 'r-', linewidth=2)
plt.plot(M, r, 'r--', linewidth=2)
plt.legend(['Якобі', 'Зейделя', 'релаксації', 'Як n0', 'Зе n0', 'ре n0'], loc='best')
plt.grid(True); plt.xlabel('$h^{-1}$'); plt.ylabel('Kit')
plt.title('Число ітерацій від h')
plt.show()
