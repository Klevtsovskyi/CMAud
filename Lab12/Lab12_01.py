"""Pr13_1: c(x,t,u)u_t = (K u')' + f(x,t,u), K=1, c=1, f=(1+u²)(1-2u); fdm_p1n1."""
import time
import numpy as np
from nonlin_heat import fdm_p1n1
from _lab12_common import a, b, u0, mn, report


def K(x, t, u):
    return 1.0


def cf(x, t, u):
    """[c(x,t,u), f(x,t,u)]"""
    return 1.0, (1 + u * u) * (1 - 2 * u)


n = 201
X = np.linspace(a, b, n)
tau, M = 2.5e-5, 16000
J = np.arange(0, M + 1, 320)
ep, sm = 1e-10, 10

t1 = time.perf_counter()
Y = fdm_p1n1(K, cf, u0, mn, X, J, tau, ep, sm)
print(f'Тривалість розрахунку ={time.perf_counter() - t1:13.3f} сек')
report(X, J, tau, Y)
