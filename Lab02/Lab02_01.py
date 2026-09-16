"""
Наближене знаходження кореня f(x) = 0.
"""

import numpy as np
import matplotlib.pyplot as plt


def m_1dihot(ffun, D, e, mki):
    """
    Знаходження кореня рівняння f(x) = 0 методом дихотомії (бісекції).

    Вхідні дані:
        ffun - функція f(x)
        D - діапазон кореня (D[0], D[1])
        e - точність (0 < e < 1)
        mki - максимальна кількість кроків ітерацій

    Вихідні дані:
        x - фінішне наближення до кореня
        k - кількість виконаних ітерацій
        err - -1 = не виконується умова D[0] < D[1] та f(D[0])*f(D[1]) <= 0
               0 = знайдено
               1 = k > mki
    """
    if e < 0 or e > 1:
        e = 1e-6

    l, r = D[0], D[1]
    fl, fr = ffun(l), ffun(r)

    if fl * fr <= 0:
        err = 1
        k = 0
        for k in range(1, mki + 1):
            x = 0.5 * (l + r)
            fx = ffun(x)
            if fl * fx <= 0:
                r, fr = x, fx
            else:
                l, fl = x, fx
            if abs(r - l) < e:
                err = 0
                break
        x = 0.5 * (l + r)
    else:
        err = -1
        k = 0
        x = None

    return x, k, err


def m_1simpit(ffun, x, e, mki):
    """
    Метод простої ітерації для розв'язання рівняння x = f(x).

    Вихідні дані:
        x, k, err - 0 = знайдено, 1 = k > mki
    """
    if e < 0 or e > 1:
        e = 1e-6

    for k in range(1, mki + 1):
        r = ffun(x)
        if abs(r - x) < e:
            return r, k, 0
        x = r

    return x, mki, 1


def m_1estef(ffun, x, e, mki):
    """
    Знаходження кореня рівняння x = f(x)
    нестаціонарним ІМ Ейткена-Стеффенсена:
        x1 = f(x_n), x2 = f(x1)
        x_{n+1} = x_n - (x_n - x1)^2 / (x_n - 2*x1 + x2)
    """
    if e < 0 or e > 1:
        e = 1e-6

    for k in range(1, mki + 1):
        x1 = ffun(x)
        x2 = ffun(x1)
        o = x - x1
        o = o * o / (x - 2.0 * x1 + x2)
        x = x - o
        if abs(o) < e:
            return x, k, 0

    return x, mki, 1


def m_1njut(ffun, x, e, mki, dffun, N):
    """
    Метод Ньютона: x_{n+1} = x_n - N * f(x_n) / f'(x_n)

    Вихідні дані:
        err - 0 = знайдено, 1 = k > mki, 2 = інші помилки
    """
    if e < 0 or e > 1:
        e = 1e-6

    if N < 1:
        return x, 0, 2

    for k in range(1, mki + 1):
        o = N * ffun(x) / dffun(x)
        x = x - o
        if abs(o) < e:
            return x, k, 0

    return x, mki, 1


def m_1secant(ffun, x, e, mki, x0):
    """
    Метод січних для знаходження розв'язку рівняння f(x) = 0.
    """
    if e < 0 or e > 1:
        e = 1e-6

    fxm = ffun(x0)
    for k in range(1, mki + 1):
        fx = ffun(x)
        o = (x - x0) / (fx - fxm) * fx
        x0 = x
        x = x - o
        if abs(o) <= e:
            return x, k, 0
        fxm = fx

    return x, mki, 1


def m_1kurch(ffun, x, e, mki, x0):
    """
    Метод М. Курчатова для знаходження розв'язку рівняння f(x) = 0.
    """
    if e < 0 or e > 1:
        e = 1e-6

    fxm = ffun(x0)
    for k in range(1, mki + 1):
        fx = ffun(x)
        o = 2 * (x - x0) / (ffun(2 * x - x0) - fxm) * fx
        x0 = x
        x = x - o
        if abs(o) <= e:
            return x, k, 0
        fxm = fx

    return x, mki, 1


tau = 0.0


def phi(x):
    # опис phi(x)
    return np.sqrt(x**3 / (x - 2))


def psi(x):
    # опис psi(x)
    return 13.1 - 1.5672 * (x - 1.2)**2


def f(x):
    # опис f(x)
    return phi(x) - psi(x)


def df(x):
    # опис f'(x)
    return 3.1344 * (x - 1.2) + (x - 3) * np.sqrt(x / (x - 2)**3)


def s(x):
    # опис s(x)
    global tau
    return x + tau * f(x)


def ds(x):
    # опис s'(x)
    global tau
    return 1.0 + tau * df(x)


def main():
    global tau

    x1 = np.linspace(-2, 4, 161)
    y1 = psi(x1)
    x2 = np.linspace(-2, 0, 81)
    y2 = phi(x2)
    x3 = np.linspace(2, 4, 41)
    x3[0] = x3[0] + 1e-2
    y3 = phi(x3)

    plt.plot(x1, y1, 'b--', linewidth=2)
    plt.plot(x2, y2, 'k', linewidth=2)
    plt.plot(x3, y3, 'k', linewidth=2)
    plt.grid(True)
    plt.legend([r'$\psi$', r'$\phi$'], loc='best')
    plt.show(block=False)

    D = [float(v) for v in input('Діапазон кореня (напр. 2 3)=').split()]
    z = float(input('Наближення до кореня='))
    plt.close('all')

    ep = 1e-8
    mki = 100  # ітераційні параметри

    X = np.linspace(D[0], D[1], 51)
    Y = df(X)
    tau = 1.0 / (np.max(np.abs(Y)) + ep)
    if df(z) > 0:
        tau = -tau

    S = ds(X)
    q = np.max(np.abs(S))

    plt.figure()
    plt.plot(x2, f(x2), 'r', linewidth=2)
    plt.plot(x3, f(x3), 'r', linewidth=2)
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('f(x)')
    plt.show()
    input('Натисніть Enter для продовження...')
    plt.close('all')

    print(f'q={q}')

    plt.figure()
    plt.plot(X, S, 'm-.', linewidth=2)
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel("s'(x)")
    plt.title("s'(x)")
    plt.show()
    input('Натисніть Enter для продовження...')
    plt.close('all')

    print(f"{'Ітераційний метод':<17s} {'X-наближен.корінь':>21s} "
          f"{'Значення f(x) в Х':>21s} {'КілІт':>5s} {'err':>3s}")

    def report(name, x, k, err):
        print(f'{name:<17s} {x:21.12e} {f(x):21.12e} {k:5d} {err:3d}')

    x, k, err = m_1dihot(f, D, ep, mki)
    report('діхотомії        ', x, k, err)

    x = z
    x, k, err = m_1simpit(s, x, ep, mki)
    report('простої ітерації ', x, k, err)

    x = z
    x, k, err = m_1estef(s, x, ep, mki)
    report('Ейтк.-Стеффенсена', x, k, err)

    x = z
    x, k, err = m_1njut(f, x, ep, mki, df, 1)
    report('Ньютона          ', x, k, err)

    x = z
    x0 = x - 0.05
    x, k, err = m_1secant(f, x, ep, mki, x0)
    report('січних           ', x, k, err)

    x = z
    x0 = x - 0.05
    x, k, err = m_1kurch(f, x, ep, mki, x0)
    report('Курчатова        ', x, k, err)


if __name__ == '__main__':
    main()
