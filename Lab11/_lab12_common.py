"""Спільні дані та виведення для Lab11_01 і Lab11_02 (модельна задача з точним розв'язком
u(x,t)=x(2-x)exp(-t))."""
import numpy as np
import matplotlib.pyplot as plt

a, b = 0, 1


def p(x, t):
    return (2 + x) * t


def f(x, t):
    return (2 + x * (x - 2) + (4 - x * x) * x * t) * np.exp(-t)


def u0(x):
    return x * (2 - x)


def u(x, t):
    return x * (2 - x) * np.exp(-t)


v = [1, 1]


def Ga(t):
    """[g(t), mu(t)] при x=a.
    КУ I роду (v[0]=0): [1, 0];  II роду (v[0]=1): [0, -2exp(-t)];  III роду: [1, -2exp(-t)]"""
    return [1, -2 * np.exp(-t)]          # III р
    # return [0, -2 * np.exp(-t)]        # II  р
    # return [1, 0]                      # I   р


def Gb(t):
    """[q(t), nu(t)] при x=b.
    КУ I роду (v[1]=0): [1, 0];  II роду: [1, exp(-t)];  III роду: [2, 2exp(-t)]"""
    return [2, 2 * np.exp(-t)]           # III р
    # return [1, np.exp(-t)]             # II  р
    # return [1, 0]                      # I   р


def report(X, J_idx, tau, U):
    J = tau * np.asarray(J_idx)
    T = J[-1]
    Mx, Mt = np.meshgrid(X, J)
    Uu = u(Mx, Mt)                         # точний розв'язок на сітці
    lim = (a, b, 0, T, -0.2, 1.2)

    def mesh(Z, title, zlabel, cmap, limits=True):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.plot_surface(Mx, Mt, Z, cmap=cmap, linewidth=0, antialiased=True)
        if limits:
            ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3]); ax.set_zlim(lim[4], lim[5])
        ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel(zlabel)
        ax.set_title(title)
        plt.show()

    mesh(U, 'U - наближений', 'U', 'cool')
    mesh(Uu, 'u(x,t) - точний', 'u', 'cool')
    E = U - Uu
    # norm(E, inf) у MATLAB для матриці - максимальна сума модулів по рядку
    print(f'Похибка ||U-u||c={np.linalg.norm(E, np.inf):13.5e}')
    print(f'max(max(|U-u|)) ={np.max(np.abs(E)):13.5e}')
    mesh(E, 'Похибка U-u(x,t)', 'U-u', 'jet', limits=False)
