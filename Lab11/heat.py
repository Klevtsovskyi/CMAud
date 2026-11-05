"""Різницеві схеми для параболічного рівняння
    u_t = u_xx - p(x,t) u + f(x,t),  x in (a,b),  p >= 0,  u(x,0) = u0(x),
з крайовими умовами
    v[0]=0: u(a,t) = mu(t)                       v[0]=1: -u_x(a,t) + g(t) u(a,t) = mu(t)
    v[1]=0: u(b,t) = nu(t)                       v[1]=1:  u_x(b,t) + q(t) u(b,t) = nu(t)
(Переклад з MATLAB: fdm_P10.m - явна схема, fdm_P1s.m - неявна схема з вагою sigm.)

Ga(t) = [g(t), mu(t)] (при v[0]=0: [1, mu(t)]), Gb(t) = [q(t), nu(t)].
Функції p(x,t), f(x,t) мають приймати масив x (векторизоване обчислення по x).
Результат U[j,k] = u(X[k], tau*J[j]).
"""
import warnings
import numpy as np
from Lab08.bvp import m_progm, m_progn


def _vec(fun, x, t):
    """fun(x,t) як масив тієї ж форми, що й x (допускається сталий результат)."""
    return np.broadcast_to(np.asarray(fun(x, t), dtype=float), np.shape(x))


def _init(u0, X, J):
    Y = np.array([u0(x) for x in X], dtype=float)
    U = np.zeros((len(J), len(X)))
    U[0] = Y
    return Y, U


def fdm_p10(p, f, u0, Ga, Gb, v, X, J, tau):
    """Явна схема (sigm = 0). Умова стійкості: 2*tau < h^2."""
    X = np.asarray(X, dtype=float)
    J = np.asarray(J, dtype=int)
    n = len(X)
    nm = n - 1
    h = X[1] - X[0]
    hh = h * h
    M = J.max()
    if 2 * tau >= hh:
        warnings.warn('FDM_P10: не виконується умова стійкості!')
    tau05 = 0.5 * tau
    gam = tau / hh
    gam2 = 2 * gam
    gamt = 2 * tau / h
    g12 = 1 - gam2
    Y, U = _init(u0, X, J)
    nJ = 1
    xin = X[1:-1]
    for j in range(1, M + 1):
        t = j * tau
        tt = t - tau05
        yold = Y
        Y = np.empty(n)
        # x = a
        if v[0] == 0:
            Y[0] = Ga(t)[1]
        else:
            gr = Ga(tt)
            pt = tau05 * p(X[0], tt)
            F = tau * f(X[0], tt) + gamt * gr[1]
            Y[0] = ((g12 - pt - gamt * gr[0]) * yold[0] + gam2 * yold[1] + F) / (1 + pt)
        # внутрішні вузли
        P = tau * _vec(p, xin, tt)
        F = tau * _vec(f, xin, tt)
        Y[1:-1] = gam * (yold[:-2] + yold[2:]) + (g12 - P) * yold[1:-1] + F
        # x = b
        if v[1] == 0:
            Y[-1] = Gb(t)[1]
        else:
            qr = Gb(tt)
            pt = tau05 * p(X[-1], tt)
            F = tau * f(X[-1], tt) + gamt * qr[1]
            Y[-1] = (gam2 * yold[nm - 1] + (g12 - pt - gamt * qr[0]) * yold[-1] + F) / (1 + pt)
        if nJ < len(J) and j == J[nJ]:
            U[nJ] = Y
            nJ += 1
    return U


def fdm_p1s(p, f, u0, Ga, Gb, v, X, J, tau, sigm):
    """Неявна схема з вагою 0 < sigm <= 1 (sigm=0.5 - Кранка-Ніколсон, 1 - чисто неявна)."""
    X = np.asarray(X, dtype=float)
    J = np.asarray(J, dtype=int)
    n = len(X)
    nm = n - 1
    h = X[1] - X[0]
    hh = h * h
    M = J.max()
    if sigm <= 0 or sigm > 1:
        sigm = 0.5
        warnings.warn('FDM_P1S: значення ваги sigm змінено на 0.5!')
    tau05 = 0.5 * tau
    s1 = 1 - sigm
    gam = tau / hh
    gam2 = 2 * gam
    gamt = 2 * tau / h
    gams = s1 * gam
    gam2s = s1 * gam2
    gam = sigm * gam
    gam2 = sigm * gam2
    g12p = 1 + gam2
    g12m = 1 - gam2s
    sgt = sigm * gamt
    s1gt = s1 * gamt
    st = sigm * tau
    s1t = s1 * tau
    Y, U = _init(u0, X, J)
    nJ = 1
    # коефіцієнти СЛАР  A(k)y(k-1) - C(k)y(k) + B(k)y(k+1) = -F(k)
    A = gam * np.ones(n)
    B = A.copy()
    C = np.ones(n)
    F = np.ones(n)
    A[0] = 0
    B[0] = 0 if v[0] == 0 else gam2
    A[-1] = 0 if v[1] == 0 else gam2
    B[-1] = 0
    xin = X[1:-1]
    for j in range(1, M + 1):
        t = j * tau
        tt = t - tau05
        if v[0] == 0:
            F[0] = Ga(t)[1]
        else:
            gr = Ga(tt)
            pt = tau05 * p(X[0], tt)
            C[0] = g12p + pt + sgt * gr[0]
            F[0] = ((g12m - pt - s1gt * gr[0]) * Y[0] + gam2s * Y[1]
                    + tau * f(X[0], tt) + gamt * gr[1])
        pt = _vec(p, xin, tt)
        C[1:-1] = g12p + st * pt
        F[1:-1] = (gams * (Y[:-2] + Y[2:]) + (g12m - s1t * pt) * Y[1:-1]
                   + tau * _vec(f, xin, tt))
        if v[1] == 0:
            F[-1] = Gb(t)[1]
        else:
            qr = Gb(tt)
            pt = tau05 * p(X[-1], tt)
            C[-1] = g12p + pt + sgt * qr[0]
            F[-1] = (gam2s * Y[nm - 1] + (g12m - pt - s1gt * qr[0]) * Y[-1]
                     + tau * f(X[-1], tt) + gamt * qr[1])
        Y, _, err = m_progm(A, C, B, F)        # монотонна прогонка
        if err:                                # в MATLAB err ігнорувалась
            Y, _ = m_progn(A, C, B, F)
        if nJ < len(J) and j == J[nJ]:
            U[nJ] = Y
            nJ += 1
    return U
