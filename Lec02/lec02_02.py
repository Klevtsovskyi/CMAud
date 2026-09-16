# -*- coding: utf-8 -*-
"""
Знаходження кореня нелінійного рівняння f(x) = 0.

Переклад з MATLAB (Le02_2.m) на Python.
Файл збережено у кодуванні utf-8.

Відповідність функцій MATLAB -> Python:
    ezplot           -> matplotlib.pyplot
    fminbnd          -> scipy.optimize.minimize_scalar(method='bounded')
    fzero(f, z)      -> scipy.optimize.newton(f, z)        (початкове наближення)
    fzero(f, [a,b])  -> scipy.optimize.brentq(f, a, b)      (відрізок з коренем)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, brentq, newton


def f(x):
    return x * x - 2 * x * np.sin(x) - 1


def F(x):
    return f(x) ** 2


def main():
    # графіки F(x) та f(x) на [-3, 3]  (аналог ezplot)
    xs = np.linspace(-3, 3, 1000)
    plt.figure()
    plt.plot(xs, F(xs), label='F(x) = f(x)^2')
    plt.plot(xs, f(xs), label='f(x)')
    plt.axhline(0, color='k', linewidth=0.5)
    plt.grid(True)
    plt.legend()
    plt.title('Le02_2: F(x) та f(x)')
    plt.show()

    d_str = input("Інтервал, де знаходиться корінь (напр.: -1 1) = ")
    D = [float(v) for v in d_str.replace(',', ' ').split()]

    z = float(input('Наближення = '))

    # Мінімізація F(x) на відрізку D
    res = minimize_scalar(F, bounds=(D[0], D[1]), method='bounded')
    x = res.x
    print('Мінімізація  : F(%.6f)=%.6e' % (x, F(x)))

    # fzero(f, z) - пошук кореня від початкового наближення z
    x = newton(f, z)
    print('Функція fzero: f(%.6f)=%.6e' % (x, f(x)))

    # fzero(f, D) - пошук кореня на відрізку [D(1), D(2)] методом Брента
    try:
        x = brentq(f, D[0], D[1], full_output=False)
        ef = 1  # код успішного завершення (аналог виводу fzero)
    except ValueError:
        # на кінцях відрізка немає зміни знаку - корінь так не знайдено
        x = float('nan')
        ef = -1
    fx = f(x) if not np.isnan(x) else float('nan')
    print('f(%.6f)=%.6e Код завершення=%2d' % (x, fx, ef))


if __name__ == '__main__':
    main()