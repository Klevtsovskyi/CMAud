import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d, CubicSpline


def view_ip(u, a, b, N, colors, title):
    nN = len(N)
    step = 5

    plt.figure(figsize=(15, 5))

    for k in range(nN):

        n = N[k]
        nm = n - 1
        Ns = nm * step + 1

        # Вузли інтерполяції
        X = np.linspace(a, b, n)

        # Згущена сітка
        x = np.linspace(a, b, Ns)

        # Значення функції у вузлах
        U = u(X)

        # ==========================================
        # 1. Поліноміальна інтерполяція
        # ==========================================

        P = np.polyfit(X, U, nm)
        Y = np.polyval(P, x)

        plt.subplot(1, 3, 1)

        plt.plot(x, Y, colors[k], label=f'N={n}')

        if k == nN - 1:
            plt.plot(X, U, 'k.', label='Вузли')
            plt.grid()
            plt.axis('equal')
            plt.xlabel('X, x')
            plt.ylabel('U, Y')
            plt.title(title)
            plt.legend()

        # ==========================================
        # 2. Лінійна інтерполяція
        # ==========================================

        linear = interp1d(
            X,
            U,
            kind='linear'
        )

        S11 = linear(x)

        plt.subplot(1, 3, 2)

        plt.plot(x, S11, colors[k], label=f'N={n}')

        if k == nN - 1:
            plt.plot(X, U, 'k.')
            plt.grid()
            plt.axis('equal')
            plt.xlabel('X, x')
            plt.ylabel('U, S11')
            plt.title(title)
            plt.legend()

        # ==========================================
        # 3. Кубічний сплайн
        # ==========================================

        spline = CubicSpline(X, U)

        S31 = spline(x)

        plt.subplot(1, 3, 3)

        plt.plot(x, S31, colors[k], label=f'N={n}')

        if k == nN - 1:
            plt.plot(X, U, 'k.')
            plt.grid()
            plt.axis('equal')
            plt.xlabel('X, x')
            plt.ylabel('U, S31')
            plt.title(title)
            plt.legend()

    plt.tight_layout()
    plt.show()


# ==========================================
# Параметри
# ==========================================

a = -1
b = 1

N = [5, 11, 17]

colors = ['b-', 'g-', 'r-']


# ==========================================
# Приклад 1
# u(t) = |t|
# ==========================================

u = lambda t: np.abs(t)

view_ip(
    u,
    a,
    b,
    N,
    colors,
    'Пр. С.Н. Бернштейна'
)


# ==========================================
# Приклад 2
# u(t) = 1 / (1 + 25t²)
# ==========================================

u = lambda t: 1 / (1 + 25 * t**2)

view_ip(
    u,
    a,
    b,
    N,
    colors,
    'Приклад К. Рунге'
)
