import numpy as np


def _adjust_n(nkf, n):
    """Корекція n: трапеції n>=2; Сімпсон n>=3 і непарне."""
    if nkf == 0:
        n = max(n, 2)
    elif nkf == 1:
        if n < 3:
            n = 3
        elif n % 2 == 0:
            n += 1
    else:
        raise ValueError('Unknown method (nkf має бути 0 або 1)')
    return n


def _grid_weights(a, b, nkf, n):
    """Вузли рівномірної сітки та коефіцієнти КФ трапецій (0) / Сімпсона (1)."""
    X = np.linspace(a, b, n)
    h = X[1] - X[0]
    P = np.empty(n)
    if nkf == 0:
        P[:] = h
        P[0] = P[-1] = 0.5 * h
    else:
        P[0] = P[-1] = h / 3
        P[1:-1:2] = 4 * h / 3
        P[2:-1:2] = 2 * h / 3
    return X, P


# ------------------------------------------------ Фредгольм 2-го роду, колокації
def fr2_mcol(lam, X, kfun, ffun):
    """y(x) = f(x) + lam * int_a^b K(x,s) y(s) ds, метод колокацій.

    X - сплайнова сітка a=X[0]<...<X[-1]=b.
    Повертає (x, Y): точки колокації (середини відрізків) і розв'язок у них.
    """
    X = np.asarray(X, dtype=float)
    if len(X) < 2:
        raise ValueError('FR2_MCOL: потрібно щонайменше 2 вузли сітки')
    N = len(X) - 1
    x = 0.5 * (X[:-1] + X[1:])
    h = np.diff(X)
    F = np.array([ffun(xi) for xi in x])
    A = np.empty((N, N))
    for i in range(N):
        for j in range(N):
            A[i, j] = -lam * h[j] * kfun(x[i], x[j])
    A += np.eye(N)
    Y = np.linalg.solve(A, F)
    return x, Y


# ------------------------------------------------ Фредгольм 2-го роду, заміна інтеграла сумою
def fr2_mk(lam, a, b, nkf, n, kfun, ffun):
    """Заміна інтеграла скінченною сумою (nkf: 0 - трапеції, 1 - Сімпсон).
    Повертає (X, Y, n)."""
    if a >= b:
        raise ValueError('FR2_MK: a >= b')
    n = _adjust_n(nkf, n)
    X, P = _grid_weights(a, b, nkf, n)
    F = np.array([ffun(x) for x in X])
    A = np.empty((n, n))
    for i in range(n):
        for k in range(n):
            A[i, k] = -lam * P[k] * kfun(X[i], X[k])
        A[i, i] += 1.0
    Y = np.linalg.solve(A, F)
    return X, Y, n


# ------------------------------------------------ Фредгольм 2-го роду, послідовні наближення
def fr2_mpn(lam, a, b, nkf, n, e, mki, kfun, ffun):
    """Метод послідовних наближень. Повертає (X, Y, n)."""
    n = _adjust_n(nkf, n)
    if a >= b:
        raise ValueError('FR2_MPN: a >= b')
    if e <= 0:
        raise ValueError('FR2_MPN: e <= 0')
    X, P = _grid_weights(a, b, nkf, n)
    FI = np.array([ffun(x) for x in X])
    Y = FI.copy()
    AK = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            AK[i, j] = lam * P[j] * kfun(X[i], X[j])

    k = 0
    F = True
    while F and k < mki:
        k += 1
        FI = AK @ FI
        Y = Y + FI
        F = bool(np.any(np.abs(FI) > e))
    return X, Y, n


# ------------------------------------------------ Фредгольм 2-го роду, вироджене ядро
def fr2_mzjd(lam, a, b, nkf, n, kfun, ffun):
    """Заміна ядра виродженим K(x,s) = sum_k A_k(x) B_k(s).

    kfun(x) повертає пару масивів (A, B) довжини m.
    Повертає (X, Y, n).
    """
    if a >= b:
        raise ValueError('FR2_MZJD: a >= b')
    n = _adjust_n(nkf, n)
    X, P = _grid_weights(a, b, nkf, n)

    cols = [kfun(x) for x in X]
    AK = np.column_stack([np.asarray(c[0], float) for c in cols])  # (m, n)
    BK = np.column_stack([np.asarray(c[1], float) for c in cols])  # (m, n)
    FF = np.array([ffun(x) for x in X])
    m = AK.shape[0]

    BP = BK * P                            # BP[i,k] = P[k]*BK[i,k]
    A = np.eye(m) - lam * (BP @ AK.T)      # A[i,j] = d_ij - lam*sum_k P_k AK[j,k] BK[i,k]
    F = BP @ FF
    C = np.linalg.solve(A, F)
    Y = FF + lam * (C @ AK)
    return X, Y, n


# ------------------------------------------------ Фредгольм-Урисон (нелінійне)
def frur_mk(a, b, nkf, n, k_dkfun, f_dffun, Y, e, mki, m):
    """int_a^b K(x,s,u(s)) ds = f(x,u(x)).

    k_dkfun(x,s,u) -> (K, dK/du);  f_dffun(x,u) -> (f, df/du).
    m: 1 - градієнта, 2 - Ньютона, 3 - спочатку градієнта (3 кроки), далі Ньютона.
    Повертає (X, Y).
    """
    if a >= b:
        raise ValueError('FRUR_MK: a >= b')
    if nkf == 0:
        if n < 2:
            raise ValueError('FRUR_MK: nkf = 0 & n < 2')
    elif nkf == 1:
        if n < 3:
            raise ValueError('FRUR_MK: nkf = 1 & n < 3')
        if n % 2 == 0:
            raise ValueError('FRUR_MK: nkf = 1 & n парне')
    else:
        raise ValueError('FRUR_MK: Unknown method')
    if m not in (1, 2, 3):
        raise ValueError('FRUR_MK: m має бути 1, 2 або 3')

    X, P = _grid_weights(a, b, nkf, n)
    Y = np.array(Y, dtype=float).ravel()

    s = 0
    Fl = True
    while Fl and s < mki:
        # матриця Якобі J і вектор F
        J = np.zeros((n, n))
        F = np.zeros(n)
        for i in range(n):
            xi = X[i]
            sk = 0.0
            for j in range(n):
                k, dk = k_dkfun(xi, X[j], Y[j])
                J[i, j] = P[j] * dk
                sk += P[j] * k
            f, df = f_dffun(xi, Y[i])
            J[i, i] -= df
            F[i] = sk - f
        if m == 1 or (m == 3 and s < 3):
            d = J.T @ F                    # метод градієнта
            v = J @ d
            d = (F @ v) / (v @ v) * d
        else:
            d = np.linalg.solve(J, F)      # метод Ньютона
        Y = Y - d
        s += 1
        Fl = np.linalg.norm(d) > e
    return X, Y


# ------------------------------------------------ Вольтерра 2-го роду
def vo2_mk(lam, a, b, n, kfun, ffun):
    """y(x) = f(x) + lam * int_a^x K(x,s) y(s) ds, КФ трапецій. Повертає (X, Y, n)."""
    if a >= b:
        raise ValueError('VO2_MK: a >= b')
    n = max(n, 2)
    X = np.linspace(a, b, n)
    h = X[1] - X[0]
    Y = np.zeros(n)
    Y[0] = ffun(X[0])          # інтеграл по [a,a] дорівнює 0
    for i in range(1, n):
        x = X[i]
        c1 = 0.0
        for k in range(i):
            p = 0.5 * h if k == 0 else h
            c1 += p * kfun(x, X[k]) * Y[k]
        c1 = ffun(x) + lam * c1
        c2 = 1 - lam * 0.5 * h * kfun(x, x)
        Y[i] = c1 / c2
    return X, Y, n
