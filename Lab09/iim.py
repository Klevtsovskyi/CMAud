"""Інтегро-інтерполяційний метод (ІІМ) для (K(x)u')' - p(x)u = -f(x).
(Переклад з MATLAB: m_iim.m, view_sol.m, view_solu.m)

Глобальні змінні MATLAB (X, n, h, h05, h25) зібрано в клас Grid.
Усі вузлові функції Kfun(j), pfun(j), ffun(j) приймають 0-based номер вузла j.
"""
import numpy as np
import matplotlib.pyplot as plt
from Lab08.bvp import m_progm, m_progn


class Grid:
    """Рівномірна сітка X на [a,b] з n вузлами."""

    def __init__(self, a, b, n):
        self.X = np.linspace(a, b, n)
        self.n = n
        self.h = self.X[1] - self.X[0]
        self.h05 = 0.5 * self.h
        self.h25 = 0.5 * self.h05

    def x_in(self, j):
        """x_j, зсунутий усередину області в крайніх вузлах (формули (4.32))."""
        x = self.X[j]
        if j == 0:
            x += self.h25
        elif j == self.n - 1:
            x -= self.h25
        return x


def m_iim(grid, Kfun, pfun, ffun, D):
    """Наближений розв'язок крайової задачі
        (K(x)u')' - p(x)u = -f(x),  K>=K0>0, p>=0, x in (a,b),
        (-D[0]*K*u' + D[1]*u)|_{x=a} = D[2],
        ( D[3]*K*u' + D[4]*u)|_{x=b} = D[5].
    D[0],D[3] = 0 і D[1],D[4] = 1 - КУ I роду; D[0],D[3] = 1 і D[1],D[4] >= 0 - КУ II, III роду.
    Повертає масив Y значень розв'язку у вузлах сітки."""
    D = list(D)
    if len(D) < 6:
        raise ValueError('M_IIM: недостатньо елементів у масиві D')
    if D[0] not in (0, 1):
        raise ValueError('M_IIM: D[0] має бути 0 або 1')
    if D[3] not in (0, 1):
        raise ValueError('M_IIM: D[3] має бути 0 або 1')
    n, h, h05 = grid.n, grid.h, grid.h05
    hh = h * h
    A = np.zeros(n)
    B = np.zeros(n)
    C = np.zeros(n)
    F = np.zeros(n)
    akp = Kfun(1)
    if D[0]:
        C[0] = akp + h * (D[1] + h05 * pfun(0))
        B[0] = akp
        F[0] = h * (D[2] + h05 * ffun(0))
    else:
        C[0] = 1
        F[0] = D[2]
    for k in range(1, n - 1):
        ak = akp
        akp = Kfun(k + 1)
        A[k], B[k] = ak, akp
        C[k] = ak + akp + hh * pfun(k)
        F[k] = hh * ffun(k)
    if D[3]:
        C[n - 1] = akp + h * (D[4] + h05 * pfun(n - 1))
        A[n - 1] = akp
        F[n - 1] = h * (D[5] + h05 * ffun(n - 1))
    else:
        C[n - 1] = 1
        F[n - 1] = D[5]
    Y, _, err = m_progm(A, C, B, F)
    if err:
        Y, _ = m_progn(A, C, B, F)
    return Y


def make_Kjump(trk, K1, K2):
    """K(x)=K1 при x<c, K2 при x>c, c = X[trk]+h/2 (trk - 0-based вузол перед розривом).
    У вузлі m=trk+1 береться напівсума."""
    m = trk + 1

    def Kf(j):
        if j < m:
            return K1
        if j == m:
            return 0.5 * (K1 + K2)
        return K2
    return Kf


def view_sol(grid, Y, txt):
    """Графік наближеного розв'язку."""
    X = grid.X
    de = 1e-1
    plt.figure()
    plt.plot(X, Y, 'b', linewidth=2)
    plt.axis([X[0] - de, X[-1] + de, Y.min() - de, Y.max() + de])
    plt.grid(True); plt.title(txt)
    plt.xlabel('x'); plt.ylabel('Y')
    plt.show()


def view_solu(grid, Y, U, ni, txt):
    """Графік Y і точного U (у ni точках), далі - абсолютна похибка."""
    X, n, h = grid.X, grid.n, grid.h
    i = np.trunc(np.linspace(1, n, ni)).astype(int) - 1
    de = 1e-1
    plt.figure()
    plt.plot(X, Y, 'b', linewidth=2)
    plt.plot(X[i], U[i], 'k*', linewidth=2)
    plt.axis([X[0] - de, X[-1] + de, U.min() - de, U.max() + de])
    plt.grid(True); plt.title(txt)
    plt.xlabel('x'); plt.ylabel('Y,u')
    plt.legend(['Y', 'u'], loc='best')
    plt.show()
    E = np.abs(Y - U)
    print(f'h^2={h * h:11.3e} Похибка в нормі C={E.max():13.5e}')
    plt.figure()
    plt.plot(X, E, 'r', linewidth=1.5)
    plt.grid(True); plt.title('Абсолютна похибка')
    plt.xlabel('x'); plt.ylabel('|Y-u(X)|')
    plt.show()
