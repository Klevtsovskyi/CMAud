"""Pr11_2: K, p неперервні, f розривна. K=1, p=0, f=-(100*sign(x)+exp(x)),
[a,b]=[-1,1], U(a)=0, U(b)=1.
Розв'язок U(x)=0.5(1-99x)+50x^2 sign(x)+exp(x)-0.5((1+x)e+(1-x)/e)."""
import numpy as np
from iim import Grid, m_iim, view_solu

a, b = -1, 1
D = [0, 1, 0, 0, 1, 1]
g = Grid(a, b, 201)
Trf = 100                        # 0-based вузол розриву f (в MATLAB: 101)


def Kf(j):
    return 1


def pf(j):
    return 0


def ff(j):
    x = g.X[j]
    if j != Trf:                 # неперервна точка, формули (4.32)
        x = g.x_in(j)
        return -100 * np.sign(x) - np.exp(x)
    # точка розриву, формула (4.33)
    return -0.5 * (np.exp(x - g.h25) + np.exp(x + g.h25))


Y = m_iim(g, Kf, pf, ff, D)


def u(x):
    return (0.5 * (1 - 99 * x) + 50 * x * x * np.sign(x) + np.exp(x)
            - 0.5 * ((1 + x) * np.exp(1) + (1 - x) * np.exp(-1)))


view_solu(g, Y, u(g.X), 7, 'K,p-неперервні,f-розривна')
