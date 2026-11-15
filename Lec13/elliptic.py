"""Ітераційні методи для різницевої задачі Діріхле для рівняння Пуассона
    -(u_xx + u_yy) = f(x,y)   (п'ятиточкова схема) на прямокутній сітці.


U[j,i] = u(X[i], Y[j])  (рядок - вузол по y, стовпець - по x).
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter


def im_jzr(m, f, muh, muv, X, Y, ro, ep, itmax):
    """Ітераційні методи: m=0 - Якобі, m=1 - Зейделя, m=2 - верхньої релаксації (параметр ro).
    f(x,y) - права частина (приймає масиви);
    muh(x) -> (u(x,Y[0]), u(x,Y[-1]));  muv(y) -> (u(X[0],y), u(X[-1],y));
    ep - точність за max|U_new - U_old|; itmax - максимальна кількість ітерацій.
    Повертає (U, it).

    Метод Зейделя й релаксації виконуються рядок за рядком: у межах рядка нове значення
    залежить від щойно знайденого лівого сусіда, тобто це лінійне рекурентне співвідношення
    y[i] = A*y[i-1] + r[i], яке розв'язується scipy.signal.lfilter (результат такий самий, як
    у послідовному циклі по i, але швидше)."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Nx, Ny = len(X), len(Y)
    Hx, Hy = X[1] - X[0], Y[1] - Y[0]
    # константи розрахункових формул
    if m == 2:
        Hx2, Hy2 = 1 / (Hx * Hx), 1 / (Hy * Hy)
        H = 2 * (Hx2 + Hy2)
        C = ro / H
        A, B, D = Hx2 * C, Hy2 * C, 1 - ro
    else:
        Hx2, Hy2 = Hx * Hx, Hy * Hy
        H = 2 * (Hx2 + Hy2)
        A, B, C, D = Hy2 / H, Hx2 / H, Hx2 * Hy2 / H, 0.0
    # граничні умови і початкове наближення
    U = np.zeros((Ny, Nx))
    c1, c2 = muh(X)
    U[0, :] = c1
    U[-1, :] = c2
    for j in range(1, Ny - 1):
        c1, c2 = muv(Y[j])
        U[j, :] = np.linspace(c1, c2, Nx)        # + лінійна інтерполяція по горизонталі
    CF = C * np.broadcast_to(np.asarray(f(X[None, :], Y[:, None]), dtype=float), (Ny, Nx))

    it, fla = 0, False
    while not fla and it < itmax:
        it += 1
        if m == 0:                               # Якобі: усі сусіди зі старого шару
            new = (A * (U[1:-1, :-2] + U[1:-1, 2:]) + B * (U[:-2, 1:-1] + U[2:, 1:-1])
                   + CF[1:-1, 1:-1])
            fla = bool(np.all(np.abs(new - U[1:-1, 1:-1]) < ep))
            U[1:-1, 1:-1] = new
        else:                                    # Зейдель / релаксація: рух знизу вгору
            fla = True
            for j in range(1, Ny - 1):
                old = U[j, 1:-1]
                r = A * U[j, 2:] + B * (U[j - 1, 1:-1] + U[j + 1, 1:-1]) + CF[j, 1:-1]
                if m == 2:
                    r = r + D * old
                new, _ = lfilter([1.0], [1.0, -A], r, zi=[A * U[j, 0]])
                if fla:
                    fla = bool(np.all(np.abs(new - old) < ep))
                U[j, 1:-1] = new
    return U, it


def _wait(s):
    """Аналог pause(s)/pause у MATLAB: s>0 - показати на s секунд, s=0 - чекати закриття вікна."""
    if s:
        plt.pause(s)
    else:
        plt.show()
    plt.close('all')


def view_u(s, X, Y, U, kc=0, UU=None):
    """Графіки різницевого розв'язку (поверхня + лінії рівня kc), а якщо задано UU - точного
    розв'язку (поверхня + лінії рівня) і похибки U-UU."""
    Mx, My = np.meshgrid(X, Y)
    plt.close('all')

    def mesh(Z, title, zlabel):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.plot_surface(Mx, My, Z, cmap='pink', linewidth=0, antialiased=True)
        ax.set_title(title); ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel(zlabel)
        _wait(s)

    def contour(Z, title):
        plt.figure()
        plt.contour(Mx, My, Z, kc, cmap='hsv')
        plt.xlabel('x'); plt.ylabel('y'); plt.title(title)
        _wait(s)

    mesh(U, 'Різницевий розв’язок', 'U')
    if kc:
        contour(U, 'Різницевий розв’язок')
    if UU is None:
        return
    mesh(UU, 'Точний розв’язок', 'u')
    if kc:
        contour(UU, 'Точний розв’язок')
    mesh(U - UU, 'Похибка U-u', 'U-u')
