"""Pr13_2: phi(u)_t = u'' + f(u), phi(u)=u, f=(1+u²)(1-2u); fdm_p1n2 (Ньютон)."""
import time
import numpy as np
from nonlin_heat import fdm_p1n2
from _lab12_common import a, b, u0, mn, report


def fi(u):
    return u


def fdf(u):
    """[f(u), df/du, dphi/du]"""
    v = 1 + u * u
    z = 2 * u
    w = 1 - z
    return v * w, z * w - 2 * v, 1.0


n = 201
X = np.linspace(a, b, n)
tau, M = 2.5e-5, 16000
J = np.arange(0, M + 1, 320)
ep, sm = 1e-10, 10

t1 = time.perf_counter()
Y = fdm_p1n2(fi, fdf, u0, mn, X, J, tau, ep, sm)
print(f'Тривалість розрахунку ={time.perf_counter() - t1:13.3f} сек')
report(X, J, tau, Y)
