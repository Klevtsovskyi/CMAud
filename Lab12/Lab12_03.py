"""Pr13_3: u_t = (K(u)u')' + f(u), K=1, f=(1+u²)(1-2u); fdm_p1n3 (предиктор-коректор)."""
import time
import numpy as np
from nonlin_heat import fdm_p1n3
from _lab12_common import a, b, u0, mn, report


def K(u):
    return 1.0


def f(u):
    return (1 + u * u) * (1 - 2 * u)


n = 201
X = np.linspace(a, b, n)
tau, M = 1e-3, 400
J = np.arange(0, M + 1, 10)

t1 = time.perf_counter()
Y = fdm_p1n3(K, f, u0, mn, X, J, tau)
print(f'Тривалість розрахунку ={time.perf_counter() - t1:13.3f} сек')
report(X, J, tau, Y)
