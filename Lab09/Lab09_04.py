"""Pr11_4: p, f неперервні, K з розривом 1-го роду.
(K u')'=0, x in (0,1); K=K1 при x<c, K2 при x>c, c=X[m]+h/2; U(0)=1, U(1)=0."""
import numpy as np
from iim import Grid, m_iim, view_solu, make_Kjump

K1, K2 = 2, 3.5
a, b = 0, 1
D = [0, 1, 1, 0, 1, 0]
g = Grid(a, b, 101)
TrK = 38                         # 0-based вузол перед розривом K (в MATLAB: 39)
c = g.X[TrK] + g.h05
t = K1 / K2
al = 1 / (t + (1 - t) * c)
be = al * t

Y = m_iim(g, make_Kjump(TrK, K1, K2), lambda j: 0, lambda j: 0, D)
U = np.where(g.X <= c, 1 - al * g.X, be * (1 - g.X))   # точний розв'язок
view_solu(g, Y, U, 21, 'p,f-неперервні,K-розривна')
