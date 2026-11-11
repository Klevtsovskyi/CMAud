"""Спільні дані Pr13_1..3: u_t = u_xx + (1+u²)(1-2u), точний розв'язок u=tan(x+t)."""
import numpy as np
import matplotlib.pyplot as plt
from nonlin_heat import view_solu_xt

a, b = 0, 1
PAUSE = 0          # тривалість показу кожного графіка, с (0 - без пауз)


def u0(x):
    return np.tan(x)


def mn(t):
    """[mu(t), nu(t)] для x=a, x=b"""
    return np.tan(t), np.tan(1 + t)


def u_exact(x, t):
    return np.tan(x + t)


def report(X, J_idx, tau, Y):
    J = tau * np.asarray(J_idx)
    U = u_exact(X[None, :], J[:, None])
    view_solu_xt(X, J, Y, U, 'Y - наближений', 'Y', 'U - точний', 'U',
                 'Похибка', 'Похибка Y-U', 'Y-U', PAUSE)
    plt.show()
