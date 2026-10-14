import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from odesolvers import rk2_ev, rk4_ev


def F(x, U):
    return np.array([np.exp(-x) - U[0]])


Y0 = [1]
X = np.linspace(0, 5, 501)
Z = np.exp(-X) * (1 + X)                 # точний розв'язок
h = X[1] - X[0]
h2 = h * h
h4 = h2 * h2


def solver(method):
    s = solve_ivp(F, [X[0], X[-1]], Y0, method=method, t_eval=X)
    return s.y[0]


def compare(Y1, Y2, n1, n2, hp, name):
    plt.plot(X, Y1, X, Y2, X, Z); plt.grid(True)
    plt.title("u'(x)=exp(-x)-u"); plt.xlabel('x'); plt.ylabel('u')
    plt.legend([n1, n2, 'точне'], loc='best')
    plt.show()
    Z1, Z2 = Z - Y1, Z - Y2
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(X, Z1); ax[0].grid(True); ax[0].set_title('u(x)-Y1')
    ax[1].plot(X, Z2); ax[1].grid(True); ax[1].set_title('u(x)-Y2')
    fig.tight_layout()
    plt.show()
    print(f"{name:>12} {'||u-Y1||c':>12} {'||u-Y2||c':>12}")
    print(f'{hp:12.4e} {np.max(np.abs(Z1)):12.4e} {np.max(np.abs(Z2)):12.4e}')


# 2-стадійний Р-К та ode23 (RK23)
compare(rk2_ev(F, X, Y0, 1)[0], solver('RK23'), 'rk2', 'ode23', h2, 'h^2')
# 4-стадійний Р-К та ode45 (RK45)
compare(rk4_ev(F, X, Y0)[0], solver('RK45'), 'rk4', 'ode45', h4, 'h^4')
