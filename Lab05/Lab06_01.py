import numpy as np
from numint import i_crectn, i_trapn, i_simpn, i_gauss


def f(x):
    return np.sin(np.pi * x)


n, a, b, ng = 7, 0, 1, 4
L = b - a
h = L / (n - 1)
I = 2 / np.pi  # точне значення int_0^1 sin(pi x) dx

Cr = i_crectn(a, b, f, n)   # центр. прямокутників
T = i_trapn(a, b, f, n)     # трапецій
S = i_simpn(a, b, f, n)     # Сімпсона
G = i_gauss(a, b, f, ng)    # Гаусса

h2 = h * h
h4 = h2 * h2
M2 = np.pi ** 2
M4 = M2 * M2
p = L * h2 * M2 / 12

print(f'{"Точне знач.=":>12}{I:19.12f}')
print(f'{"Назва КФ":>12} {"Наближене значення":>12} '
      f'{"Абсолютна  похибка":>12} {"Теор.оцінк.похибки":>12}')


def row(name, val, err, est):
    print(f'{name:>12}{val:19.12f}{err:19.10e}{est:19.10e}')


row(f'центр.пр n={n}', Cr, abs(I - Cr), p / 2)
row(f'трапецій n={n}', T, abs(I - T), p)
row(f'Сімпсона n={n}', S, abs(I - S), L * h4 * M4 / 180)
row(f'Гаусса   n={ng}', G, abs(I - G), 1e-8 * M4 * M4)
