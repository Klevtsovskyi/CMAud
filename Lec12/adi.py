"""Двовимірне параболічне рівняння u_t = u_x1x1 + u_x2x2 + f(x1,x2,t) з КУ I роду:
метод змінних напрямків (схема Пісмена-Рекфорда, ППС).
(Переклад з MATLAB: fmd_P2PR.m, view_solu.m [2D-версія])

Масиви розв'язку мають розмір (N2, N1): Y[k2,k1] = u(X1[k1], X2[k2]).
Функції f, muh, muv, u0 мають приймати масиви x1,x2 (векторизоване обчислення).
"""
import numpy as np
import matplotlib.pyplot as plt


def m_progm_multi(a, c, b, F):
    """Монотонна прогонка для багатьох правих частин одразу.
    Коефіцієнти a,c,b (довжини n) спільні, F має розмір (n, m): кожен стовпець - своя СЛАР
        a[k]*y[k-1] - c[k]*y[k] + b[k]*y[k+1] = -F[k].
    Повертає масив y розміру (n, m). Результат збігається з m_progm для кожного стовпця."""
    a, c, b, F = (np.asarray(v, dtype=float) for v in (a, c, b, F))
    n = len(c)
    if np.any(np.abs(c) < np.abs(a) + np.abs(b)):
        raise ValueError('m_progm_multi: не виконано умов монотонної прогонки')
    alfa = np.zeros(n)
    beta = np.zeros_like(F)
    for k in range(n - 1):
        z = c[k] - a[k] * alfa[k]
        alfa[k + 1] = b[k] / z
        beta[k + 1] = (F[k] + a[k] * beta[k]) / z
    y = np.zeros_like(F)
    y[n - 1] = (F[n - 1] + a[n - 1] * beta[n - 1]) / (c[n - 1] - a[n - 1] * alfa[n - 1])
    for k in range(n - 2, -1, -1):
        y[k] = alfa[k + 1] * y[k + 1] + beta[k + 1]
    return y


def fmd_p2pr(f, muh, muv, u0, X1, X2, tau, M, J=()):
    """Лінійна мішана задача для 2D параболічного рівняння з КУ I роду (ППС).
    Вхідні дані:
      f(x1,x2,t); u0(x1,x2);
      muh(x1,t) -> (mu(x1,X2[0],t), mu(x1,X2[-1],t))   - "горизонталі";
      muv(x2,t) -> (mu(X1[0],x2,t), mu(X1[-1],x2,t))   - "вертикалі";
      X1, X2 - рівномірні сітки; tau; M (T=M*tau);
      J - зростаючий масив номерів шарів, розв'язок на яких записується у файли P2PR_<j>.npy
          (порожній - нічого не записується).
    Повертає Y (N2,N1) на шарі j=M.
    """
    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)
    N1, N2 = len(X1), len(X2)
    h1, h2 = X1[1] - X1[0], X2[1] - X2[0]
    tau05 = 0.5 * tau
    ga1 = tau05 / (h1 * h1)
    s1 = 1 - 2 * ga1
    c1 = 1 + 2 * ga1
    ga2 = tau05 / (h2 * h2)
    s2 = 1 - 2 * ga2
    c2 = 1 + 2 * ga2
    # коефіцієнти СЛАР, що не залежать від j (A y[k-1] - C y[k] + B y[k+1] = -F)
    A1 = np.r_[0.0, ga1 * np.ones(N1 - 2), 0.0]
    B1 = A1.copy()
    C1 = np.ones(N1)
    C1[1:-1] = c1
    A2 = np.r_[0.0, ga2 * np.ones(N2 - 2), 0.0]
    B2 = A2.copy()
    C2 = np.ones(N2)
    C2[1:-1] = c2
    # початкові умови
    G1, G2 = np.meshgrid(X1, X2)                  # форма (N2,N1)
    Y = np.asarray(u0(G1, G2), dtype=float).copy()
    Y1 = np.zeros((N2, N1))
    J = list(J)
    jout = 0
    if J and J[0] == 0:
        np.save('P2PR_0.npy', Y)
        jout = 1
    x1i, x2i = X1[1:-1], X2[1:-1]
    for j in range(1, M + 1):
        t = j * tau
        tt = t - tau05
        tm = (j - 1) * tau
        # --- розв'язок на шарі j-1/2 (неявно по x1) ---
        ga_, gb_ = muh(X1, tt)
        Y1[0, :] = ga_
        Y1[-1, :] = gb_
        Va, Vb = muv(X2, t)                       # muv на шарі t
        Wa, Wb = muv(X2, tm)                      # muv на шарі t-1
        F1 = np.empty((N1, N2 - 2))               # стовпець - розрахунковий рядок k2
        # крайові умови проміжного шару при k1=1, N1
        for row, (W, V) in ((0, (Wa, Va)), (N1 - 1, (Wb, Vb))):
            F1[row, :] = 0.5 * ((W[1:-1] + V[1:-1]) + ga2 * (
                (W[:-2] - V[:-2]) - 2 * (W[1:-1] - V[1:-1]) + (W[2:] - V[2:])))
        fg = np.asarray(f(x1i[None, :], x2i[:, None], tt), dtype=float) * np.ones((N2 - 2, N1 - 2))
        G = ga2 * (Y[:-2, 1:-1] + Y[2:, 1:-1]) + s2 * Y[1:-1, 1:-1] + tau05 * fg
        F1[1:-1, :] = G.T
        Y1[1:-1, :] = m_progm_multi(A1, C1, B1, F1).T
        # --- розв'язок на шарі j (неявно по x2) ---
        Y[:, 0] = Va
        Y[:, -1] = Vb
        ca, cb = muh(x1i, t)
        F2 = np.empty((N2, N1 - 2))               # стовпець - розрахунковий стовпець k1
        F2[0, :] = ca
        F2[-1, :] = cb
        fg = np.asarray(f(x1i[None, :], x2i[:, None], tt), dtype=float) * np.ones((N2 - 2, N1 - 2))
        F2[1:-1, :] = ga1 * (Y1[1:-1, :-2] + Y1[1:-1, 2:]) + s1 * Y1[1:-1, 1:-1] + tau05 * fg
        Y[:, 1:-1] = m_progm_multi(A2, C2, B2, F2)
        if jout < len(J) and J[jout] == j:
            np.save(f'P2PR_{j}.npy', Y)
            jout += 1
    return Y


def view_solu2d(X1, X2, Y, z1, ti1, pnt=0, s=0, U=None, z2='', ti2='',
                txtp='', z3='', ti3=''):
    """Графіки: наближений розв'язок Y; якщо задано U - ще точний розв'язок і похибка Y-U.
    pnt - друкувати масив Y; параметр s (пауза) збережено для сумісності."""
    N1, N2 = len(X1), len(X2)
    mY, MY = np.min(Y), np.max(Y)
    if pnt:
        for k2 in range(N2 - 1, -1, -1):
            print(''.join(f'{v:11.2e}' for v in Y[k2, :]))
    Mx1, Mx2 = np.meshgrid(X1, X2)
    lim = (X1[0], X1[-1], X2[0], X2[-1], mY - 0.1, MY + 0.1)

    def mesh(Z, title, zlabel, cmap, limits=True):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.plot_surface(Mx1, Mx2, Z, cmap=cmap, linewidth=0, antialiased=True)
        if limits:
            ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3]); ax.set_zlim(lim[4], lim[5])
        ax.set_xlabel('x1'); ax.set_ylabel('x2'); ax.set_zlabel(zlabel)
        ax.set_title(title)
        plt.show()

    mesh(Y, ti1, z1, 'cool')
    if U is None:
        return
    mesh(U, ti2, z2, 'cool')
    YU = Y - U
    # norm(YU, inf) у MATLAB для матриці - максимальна сума модулів по рядку
    print(f'{txtp} ||Y-U||c={np.linalg.norm(YU, np.inf):13.5e}')
    print(f'max(max(|Y-U|)) ={np.max(np.abs(YU)):13.5e}')
    mesh(YU, ti3, z3, 'jet', limits=False)
