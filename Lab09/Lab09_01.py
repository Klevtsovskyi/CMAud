"""Pr11_1: K, p, f неперервні. K=1+x, p=exp(-x), f=sin(x)(1+x+exp(-x))-cos(x),
[a,b]=[1,3], U(a)=sin(a), U(b)=sin(b); розв'язок U(x)=sin(x)."""
import numpy as np
from iim import Grid, m_iim, view_solu

a, b = 1, 3
sa, sb = np.sin(a), np.sin(b)
# Варіанти крайових умов (в оригіналі діє останній - III р., III р.):
# [0,1,sa,0,1,sb]                               I р.,   I р.
# [1,0,-(1+a)*cos(a),0,1,sb]                    II р.,  I р.
# [1,2,-(1+a)*cos(a)+2*sa,0,1,sb]               III р., I р.
# [0,1,sa,1,0,(1+b)*cos(b)]                     I р.,   II р.
# [1,0,-(1+a)*cos(a),1,0,(1+b)*cos(b)]          II р.,  II р.
# [1,2,-(1+a)*cos(a)+2*sa,1,0,(1+b)*cos(b)]     III р., II р.
# [0,1,sa,1,3,(1+b)*cos(b)+3*sb]                I р.,   III р.
# [1,0,-(1+a)*cos(a),1,3,(1+b)*cos(b)+3*sb]     II р.,  III р.
D = [1, 2, -(1 + a) * np.cos(a) + 2 * sa, 1, 3, (1 + b) * np.cos(b) + 3 * sb]

g = Grid(a, b, 201)


def Kf(j):                       # K(x)=1+x у півцілому вузлі x_j - h/2
    return 1 + (g.X[j] - g.h05)


def pf(j):                       # p(x)=exp(-x)
    return np.exp(-g.x_in(j))


def ff(j):                       # f(x)
    x = g.x_in(j)
    return np.sin(x) * (1 + x + np.exp(-x)) - np.cos(x)


Y = m_iim(g, Kf, pf, ff, D)
U = np.sin(g.X)                  # точний розв'язок
view_solu(g, Y, U, 7, 'K,p,f-неперервні')
