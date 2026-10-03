import numpy as np
from scipy.integrate import quad
from numint import i_gauss, i_cheb, i_markov

a = 2 / 23
a1 = 1 - a * a


def f(x):
    return np.sqrt(1 - a1 * np.sin(x) ** 2)


tol = 1e-10
A, B = 0, 0.5 * np.pi

N = [2, 4, 6, 7]
E = 4 * quad(f, A, B, epsabs=tol, epsrel=tol)[0]
print(f'E= {E:.11e} - вбудована процедура (scipy.integrate.quad)')
print(' n КФ Гаусса\t    КФ Чебишева\t     КФ Маркова')
for n in N:
    g = 4 * i_gauss(A, B, f, n)
    c = 4 * i_cheb(A, B, f, n)
    m = 4 * i_markov(A, B, f, n)
    print(f'{n:2d} {g:.9e} {c:.9e} {m:.9e}')
