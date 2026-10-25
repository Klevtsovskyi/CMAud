"""Pr11_3: K, f неперервні, p розривна.
u''-p(x)u=-sin(pi x), x in (0,1); p=p1 при x<c, p2 при x>c, c=X[m];
U(0)=1, U'(1)+3U(1)=4."""
import numpy as np
from iim import Grid, m_iim, view_sol

p1, p2 = 2, 10
a, b = 0, 1
D = [0, 1, 1, 1, 3, 4]
g = Grid(a, b, 101)
Trp = 30                         # 0-based вузол розриву p (в MATLAB: 31)


def Kf(j):
    return 1


def pf(j):
    if j == Trp:                 # точка розриву, формула (4.33)
        return 0.5 * (p1 + p2)
    return p1 if j < Trp else p2


def ff(j):
    return np.sin(np.pi * g.x_in(j))


Y = m_iim(g, Kf, pf, ff, D)
view_sol(g, Y, 'K,f-неперервні,p-розривна')
