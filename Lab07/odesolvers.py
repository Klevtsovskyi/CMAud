"""Методи розв'язання задачі Коші U'(x)=F(x,U), U(X[0])=U0 на фіксованій сітці X.

Повертають масив Y розміру (dim, len(X)): k-й стовпець - розв'язок у X[k].
(Переклад з MATLAB: Euler_ev.m, Euler_nev.m, symm_nev.m, rk2_ev.m, rk4_ev.m)
"""
import numpy as np
from scipy.optimize import fsolve


def _init(X, U0):
    U0 = np.atleast_1d(np.asarray(U0, dtype=float)).ravel()
    Y = np.zeros((U0.size, len(X)))
    Y[:, 0] = U0
    return Y


def _F(F, x, z):
    return np.atleast_1d(np.asarray(F(x, z), dtype=float)).ravel()


def euler_ev(F, X, U0):
    """Явний метод Ейлера."""
    Y = _init(X, U0)
    for k in range(len(X) - 1):
        x = X[k]
        h = X[k + 1] - x
        z = Y[:, k]
        Y[:, k + 1] = z + h * _F(F, x, z)
    return Y


def _solve(fun, y0, epsx, epsf, k, x):
    sol, info, ier, msg = fsolve(fun, y0, xtol=epsx, full_output=True)
    # при дуже малому xtol scipy може не "підтвердити" збіжність (ier=5),
    # тому додатково перевіряємо нев'язку
    if ier != 1 and np.max(np.abs(fun(sol))) > max(epsf, 1e-10):
        raise RuntimeError(f'fsolve не збігся на кроці k={k}, x={x}: {msg}')
    return sol


def euler_nev(F, X, U0, epsx=1e-12, epsf=1e-12):
    """Неявний метод Ейлера (нелінійне рівняння на кожному кроці - fsolve).
    epsf - допуск на нев'язку нелінійного рівняння."""
    Y = _init(X, U0)
    for k in range(len(X) - 1):
        x = X[k + 1]
        h = x - X[k]
        z = Y[:, k]
        Y[:, k + 1] = _solve(lambda T: T - z - h * _F(F, x, T), z, epsx, epsf, k + 1, x)
    return Y


def symm_nev(F, X, U0, epsx=1e-12, epsf=1e-12):
    """Симетрична (трапецій) неявна схема."""
    Y = _init(X, U0)
    for k in range(len(X) - 1):
        x = X[k + 1]
        h = 0.5 * (x - X[k])
        z = Y[:, k]
        zf = z + h * _F(F, X[k], z)
        Y[:, k + 1] = _solve(lambda T: T - zf - h * _F(F, x, T), z, epsx, epsf, k + 1, x)
    return Y


def rk2_ev(F, X, U0, gam=1.0):
    """Явний 2-стадійний метод Рунге-Кутта 2-го порядку; 0 < gam <= 1
    (gam=1 - предиктор-коректор / середніх точок, gam=0.5 - Гюна)."""
    if gam <= 0 or gam > 1:
        gam = 1.0
    Y = _init(X, U0)
    A1, A2, g2 = 1 - gam, gam, 0.5 / gam
    for k in range(len(X) - 1):
        x = X[k]
        h = X[k + 1] - x
        hg = h * g2
        z = Y[:, k]
        K1 = _F(F, x, z)
        K2 = _F(F, x + hg, z + hg * K1)
        Y[:, k + 1] = z + h * (A1 * K1 + A2 * K2)
    return Y


def rk4_ev(F, X, U0):
    """Класичний явний 4-стадійний метод Рунге-Кутта 4-го порядку."""
    Y = _init(X, U0)
    a1, a2 = 1 / 6, 1 / 3
    for k in range(len(X) - 1):
        x, xp = X[k], X[k + 1]
        h = xp - x
        g = 0.5 * h
        z = Y[:, k]
        K1 = _F(F, x, z)
        xt = x + g
        K2 = _F(F, xt, z + g * K1)
        K3 = _F(F, xt, z + g * K2)
        K4 = _F(F, xp, z + h * K3)
        Y[:, k + 1] = z + h * (a1 * K1 + a2 * K2 + a2 * K3 + a1 * K4)
    return Y
