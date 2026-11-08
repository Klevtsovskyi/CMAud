"""Le13_1: метод змінних напрямків (ППС) для модельної 2D задачі
u_t = Δu + f,  f = exp(-t)(4+x1²+x2²) - 7,  u(x,0)=x1²+x2²,
точний розв'язок u = t + (x1²+x2²)(2-exp(-t)),  x1 in [0,1], x2 in [0,2]."""
import time
import numpy as np
from adi import fmd_p2pr, view_solu2d

a1, b1, a2, b2 = 0, 1, 0, 2


def f(x1, x2, t):                  # права частина
    return np.exp(-t) * (4 + x1 * x1 + x2 * x2) - 7


def u0(x1, x2):                    # початкова умова
    return x1 * x1 + x2 * x2


def muh(x1, t):
    """[mu(x1,a2,t), mu(x1,b2,t)]"""
    xx = x1 * x1
    e = 2 - np.exp(-t)
    return t + xx * e, t + (xx + 4) * e


def muv(x2, t):
    """[mu(a1,x2,t), mu(b1,x2,t)]"""
    xx = x2 * x2
    e = 2 - np.exp(-t)
    return t + xx * e, t + (1 + xx) * e


N1, N2 = 21, 41
X1 = np.linspace(a1, b1, N1)
X2 = np.linspace(a2, b2, N2)
tau, M = 0.0025, 2000
J = []

t1 = time.perf_counter()
Y = fmd_p2pr(f, muh, muv, u0, X1, X2, tau, M, J)
print(f'Тривалість розрахунку ={time.perf_counter() - t1:13.3f} сек')

T = tau * M
G1, G2 = np.meshgrid(X1, X2)
U = T + (G1 ** 2 + G2 ** 2) * (2 - np.exp(-T))       # точний розв'язок
view_solu2d(X1, X2, Y, 'Y', 'Y - наближений', 0, 5, U, 'U', 'U - точний',
            'Похибка', 'Y-U', 'Похибка Y-U')
