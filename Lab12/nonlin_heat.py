"""Різницеві схеми для нелінійних параболічних рівнянь (одновимірних) з КУ I роду
    u(a,t) = mu(t),  u(b,t) = nu(t),  u(x,0) = u0(x).

  fdm_p1n1:  c(x,t,u) u_t = (K(x,t,u) u_x)_x + f(x,t,u)   (прості ітерації на кожному шарі)
  fdm_p1n2:  (phi(u))_t = u_xx + f(u)                     (метод Ньютона на кожному шарі)
  fdm_p1n3:  u_t = (K(u) u_x)_x + f(u)                    (двокрокова схема предиктор-коректор)

Результат U[j,k] = u(X[k], tau*J[j]). Усі функції K, cf, fdf, fi, f мають приймати масиви
(векторизоване обчислення по x); сталі значення допускаються.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded


def _b(v, ref):
    """v як масив тієї ж форми, що й ref (допускається скаляр)."""
    return np.broadcast_to(np.asarray(v, dtype=float), np.shape(ref))


def _setup(u0, X, J):
    X = np.asarray(X, dtype=float)
    J = np.asarray(J, dtype=int)
    Y = np.array([u0(x) for x in X], dtype=float)
    U = np.zeros((len(J), len(X)))
    U[0] = Y
    return X, J, Y, U


def _tri(A, C, B, F):
    """Розв'язок СЛАР  A[k]y[k-1] - C[k]y[k] + B[k]y[k+1] = -F[k]  (LAPACK, solve_banded).
    Математично збігається з m_progm / m_progn, але у ~10 разів швидше в циклах по часу."""
    ab = np.zeros((3, len(C)))
    ab[0, 1:] = B[:-1]
    ab[1, :] = -C
    ab[2, :-1] = A[1:]
    return solve_banded((1, 1), ab, -np.asarray(F, dtype=float))


def fdm_p1n1(K, cf, u0, mn, X, J, tau, ep, smax):
    """c(x,t,u)u_t=(K(x,t,u)u')'+f(x,t,u). K(x,t,u) -> K; cf(x,t,u) -> (c, f);
    mn(t) -> (mu(t), nu(t)); ep - точність ітерацій по s; smax - макс. кількість ітерацій."""
    X, J, Y, U = _setup(u0, X, J)
    n = len(X)
    h = X[1] - X[0]
    gam05 = 0.5 * tau / (h * h)
    M = J.max()
    nJ = 1
    A = np.zeros(n)
    B = np.zeros(n)
    F = np.zeros(n)
    C = np.ones(n)
    xin = X[1:-1]
    for j in range(1, M + 1):
        t = j * tau
        Yp = Y
        F[0], F[-1] = mn(t)
        for s in range(smax):
            Ys = Y
            Kv = _b(K(X, t, Y), Y)
            c_, f_ = cf(xin, t, Y[1:-1])
            A[1:-1] = gam05 * (Kv[:-2] + Kv[1:-1])
            B[1:-1] = gam05 * (Kv[1:-1] + Kv[2:])
            C[1:-1] = A[1:-1] + B[1:-1] + c_
            F[1:-1] = c_ * Yp[1:-1] + tau * f_
            Y = _tri(A, C, B, F)
            if np.max(np.abs(Y - Ys)) <= ep:
                break
        if nJ < len(J) and j == J[nJ]:
            U[nJ] = Y
            nJ += 1
    return U


def fdm_p1n2(fi, fdf, u0, mn, X, J, tau, ep, smax):
    """phi(u)_t = u'' + f(u). fi(u) = phi(u); fdf(u) -> (f(u), df/du, dphi/du);
    ітерації Ньютона на кожному шарі."""
    X, J, Y, U = _setup(u0, X, J)
    n = len(X)
    h = X[1] - X[0]
    gam = tau / (h * h)
    gam2 = 2 * gam
    M = J.max()
    nJ = 1
    A = gam * np.ones(n)
    A[0] = A[-1] = 0
    B = A.copy()
    C = np.ones(n)
    F = np.ones(n)
    for j in range(1, M + 1):
        t = j * tau
        Yp = Y
        mu, nu = mn(t)
        for s in range(smax):
            Ys = Y
            F[0] = mu - Ys[0]
            F[-1] = nu - Ys[-1]
            y = Ys[1:-1]
            fv, dfv, dphi = fdf(y)
            C[1:-1] = gam2 - tau * dfv + dphi
            F[1:-1] = (gam * (Ys[:-2] - 2 * y + Ys[2:]) + tau * fv
                       + _b(fi(Yp[1:-1]), y) - _b(fi(y), y))
            dY = _tri(A, C, B, F)
            x = np.max(np.abs(dY))
            Y = Ys + dY
            if x <= ep:
                break
        if nJ < len(J) and j == J[nJ]:
            U[nJ] = Y
            nJ += 1
    return U


def fdm_p1n3(K, f, u0, mn, X, J, tau):
    """u_t = (K(u)u')' + f(u): схема з двома напівкроками (по t-tau/2 і t)."""
    X, J, Y, U = _setup(u0, X, J)
    n = len(X)
    h = X[1] - X[0]
    gam05 = 0.5 * tau / (h * h)
    tau05 = 0.5 * tau
    tau2 = 2 * tau
    M = J.max()
    nJ = 1
    A = np.zeros(n)
    B = np.zeros(n)
    F = np.zeros(n)
    C = np.ones(n)
    for j in range(1, M + 1):
        # СЛАР на проміжному шарі
        t = (j - 0.5) * tau
        F[0], F[-1] = mn(t)
        Kv = _b(K(Y), Y)
        A[1:-1] = gam05 * (Kv[:-2] + Kv[1:-1])
        B[1:-1] = gam05 * (Kv[1:-1] + Kv[2:])
        C[1:-1] = A[1:-1] + B[1:-1] + 2
        F[1:-1] = 2 * Y[1:-1] + tau * _b(f(Y[1:-1]), Y[1:-1])
        Yt = _tri(A, C, B, F)
        # СЛАР на шарі j
        t = t + tau05
        F[0], F[-1] = mn(t)
        Kv = _b(K(Yt), Yt)
        ai = gam05 * (Kv[:-2] + Kv[1:-1])
        bi = gam05 * (Kv[1:-1] + Kv[2:])
        A[1:-1], B[1:-1] = ai, bi
        x = ai + bi
        C[1:-1] = 2 + x
        F[1:-1] = (ai * Y[:-2] + (2 - x) * Y[1:-1] + bi * Y[2:]
                   + tau2 * _b(f(Yt[1:-1]), x))
        Y = _tri(A, C, B, F)
        if nJ < len(J) and j == J[nJ]:
            U[nJ] = Y
            nJ += 1
    return U


def view_solu_xt(X, J, Y, U, ti1, z1, ti2, z2, txtp, ti3, z3, s=0):
    """Три поверхні (x,t): наближений розв'язок Y, точний U і похибка Y-U.
    Паузи s секунд між графіками (s=0 - без пауз); вікна лишаються відкритими -
    наприкінці викличте plt.show()."""
    X = np.asarray(X)
    J = np.asarray(J)
    Mx, Mt = np.meshgrid(X, J)
    mU, MU = U.min(), U.max()
    lim = (X[0], X[-1], J[0], J[-1], mU - 0.1, MU + 0.1)

    def mesh(Z, title, zlabel, cmap, limits=True, view=None):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.plot_surface(Mx, Mt, Z, cmap=cmap, linewidth=0, antialiased=True)
        if limits:
            ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3]); ax.set_zlim(lim[4], lim[5])
        if view:
            ax.view_init(elev=view[1], azim=view[0])
        ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_zlabel(zlabel); ax.set_title(title)
        if s:
            plt.pause(s)

    plt.close('all')
    mesh(Y, ti1, z1, 'cool')
    mesh(U, ti2, z2, 'cool')
    YU = Y - U

    print(f'{txtp} ||Y-U||c={np.linalg.norm(YU, np.inf):13.5e}')
    print(f'max(max(|Y-U|)) ={np.max(np.abs(YU)):13.5e}')
    mesh(YU, ti3, z3, 'jet', limits=False, view=(-145, 35))
