"""Lab13_01: задача Діріхле/Неймана для рівняння Пуассона  -(u_xx+u_yy) = f  в області,
складеній з прямокутників і трикутника (метод верхньої релаксації):
   область: {0<=x<=1, 0<=y<=1} без квадрата-вирізу {0<=x<0.3, 0<=y<0.3} і над прямою y=2-2x.
   D1: y=0.3, 0<=x<=0.3;  D3: y=0, 0.3<=x<=1;  D5: y=1, 0<=x<=0.5;  D4: y=2-2x, 0.5<=x<=1
       - КУ I роду;  D2: x=0.3, 0<y<0.3  та  D6: x=0, 0.3<y<1 - КУ II роду (u_x задано).
Точний розв'язок u = 15 x(1-x) y(1-y).
Масив U[j,i] = u(X[i], Y[j])."""
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from scipy.signal import lfilter

a = b = 1


def f(x, y):
    return 30 * (x * (1 - x) + y * (1 - y))


def mux11(x):                      # D1 = {0<=x<=0.3, y=0.3}
    return 3.15 * x * (1 - x)


def mux12(x):                      # D3 = {0.3<=x<=1, y=0}
    return 0.0 * x


def mux21(x):                      # D5 = {0<=x<=0.5, y=1}
    return 0.0 * x


def muxy(x):                       # D4 = {0.5<=x<=1, y=2-2x}
    xx = 2 * x
    return 15 * xx * (1 - x) ** 2 * (xx - 1)


def muy11(y):                      # D6 = {0.3<y<1, x=0}: u_x
    return 15 * y * (1 - y)


def muy12(y):                      # D2 = {0<y<0.3, x=0.3}: u_x
    return 6 * y * (1 - y)


def solve(Nx=201, Ny=101, Nx11=61, Nx21=101, Ny11=31, ro=1.965, ep=1e-8, itmax=30000):
    """Метод верхньої релаксації (Зейдель з параметром ro, послідовний обхід знизу вгору
    і зліва направо, як в оригіналі). Номери Nx11, Nx21, Ny11 - 1-based, як у MATLAB.
    Повертає (X, Y, U, it)."""
    X = np.linspace(0, a, Nx)
    Y = np.linspace(0, b, Ny)
    Hx, Hy = X[1] - X[0], Y[1] - Y[0]
    HxH, HyH = 1 / (Hx * Hx), 1 / (Hy * Hy)
    C = ro * 0.5 / (HxH + HyH)
    A, B, D = HxH * C, HyH * C, 1 - ro
    A2 = 2 * A
    # крайові умови та початкове наближення
    U = np.zeros((Ny, Nx))
    U[Ny11 - 1, :Nx11] = mux11(X[:Nx11])
    U[0, Nx11 - 1:] = mux12(X[Nx11 - 1:])
    U[Ny - 1, :Nx21] = mux21(X[:Nx21])
    for j1 in range(2, Ny):                      # D4: i = Nx-j+1 (1-based)
        U[j1 - 1, Nx - j1] = muxy(X[Nx - j1])
    CF = C * f(X[None, :], Y[:, None])
    # опис рядків, що обчислюються: (j, початок, кінець (включно), особливий перший вузол, u_x)
    rows = []
    for j1 in range(2, Ny11):                    # x=0.3 & 0<y<0.3 та 0.3<x<1
        rows.append((j1 - 1, Nx11 - 1, Nx - j1 - 1, True, muy12(Y[j1 - 1])))
    rows.append((Ny11 - 1, Nx11, Nx - Ny11 - 1, False, 0.0))       # y=0.3, 0.3<x<1
    for j1 in range(Ny11 + 1, Ny):               # x=0 & 0.3<y<1 та 0<x<1
        rows.append((j1 - 1, 0, Nx - j1 - 1, True, muy11(Y[j1 - 1])))

    it, fla = 0, False
    while not fla and it < itmax:
        fla = True
        it += 1
        for j, s, e, special, mu in rows:
            sl = slice(s, e + 1)
            old = U[j, sl].copy()
            r = (A * U[j, s + 1:e + 2] + B * (U[j - 1, sl] + U[j + 1, sl])
                 + D * old + CF[j, sl])
            if special:      # вузол з умовою II роду: фіктивний лівий вузол через u_x
                r[0] = (A2 * (U[j, s + 1] - Hx * mu) + B * (U[j - 1, s] + U[j + 1, s])
                        + D * old[0] + CF[j, s])
                zi = [0.0]
            else:
                zi = [A * U[j, s - 1]]
            # у межах рядка нове значення залежить від нового лівого сусіда:
            # y[i] = A*y[i-1] + r[i] - лінійна рекурентність
            new, _ = lfilter([1.0], [1.0, -A], r, zi=zi)
            if fla:
                fla = bool(np.all(np.abs(new - old) < ep))
            U[j, sl] = new
    return X, Y, U, it


def exact_on_domain(X, Y, Nx, Ny, Nx11, Ny11):
    """Точний розв'язок в області (поза областю - 0)."""
    I1, J1 = np.meshgrid(np.arange(1, Nx + 1), np.arange(1, Ny + 1))
    mask = ((J1 >= Ny11) | (I1 >= Nx11)) & (I1 <= Nx - J1 + 1)
    Xg, Yg = np.meshgrid(X, Y)
    return np.where(mask, 15 * Xg * (1 - Xg) * Yg * (1 - Yg), 0.0)


def main():
    Nx, Ny, Nx11, Nx21, Ny11 = 201, 101, 61, 101, 31
    t0 = time.perf_counter()
    X, Y, U, it = solve(Nx, Ny, Nx11, Nx21, Ny11, ro=1.965, ep=1e-8, itmax=30000)
    print(f'Виконано ітерацій ={it:5d}')
    print(f'Тривалість розрах.={time.perf_counter() - t0:9.3f} сек')

    UU = exact_on_domain(X, Y, Nx, Ny, Nx11, Ny11)
    Mx, My = np.meshgrid(X, Y)
    d = 5e-4
    lim = (-d, a + d, -d, b + d, -d, 1 + d)
    ls = LightSource(azdeg=315, altdeg=45)

    def surf(Z, title, zlabel):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        rgb = ls.shade(Z, cmap=plt.cm.pink, vert_exag=0.1, blend_mode='soft')
        ax.plot_surface(Mx, My, Z, facecolors=rgb, linewidth=0, antialiased=False, shade=False)
        # ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3]); ax.set_zlim(lim[4], lim[5])
        ax.view_init(elev=20, azim=40 - 22.5)
        ax.set_title(title); ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel(zlabel)
        plt.show()

    def contour(Z, title):
        plt.figure()
        plt.contour(Mx, My, Z, 15, cmap='hsv')
        plt.xlabel('x'); plt.ylabel('y'); plt.title(title)
        plt.show()

    surf(U, 'Різницевий розв’язок', 'U')
    contour(U, 'Різницевий розв’язок')
    surf(UU, 'Точний розв’язок', 'u')
    contour(UU, 'Точний розв’язок')
    E = U - UU
    # norm(E, inf) у MATLAB для матриці - максимальна сума модулів по рядку
    print(f'Похибка U-u  ={np.linalg.norm(E, np.inf):13.5e}')
    print(f'max(max|U-u|)={np.max(np.abs(E)):13.5e}')
    surf(E, 'Похибка U-u', 'U-u')


if __name__ == '__main__':
    main()
