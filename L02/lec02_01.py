# -*- coding: utf-8 -*-
"""
Контрольний розрахунок для системи нелінійних рівнянь
    F(X) = 0,   де X_точне = 1.

Точний переклад з MATLAB на Python:
    ns_main.m, m_ns_njut.m, m_ns_njutg.m, m_ns_broy.m, m_ns_broym.m,
    m_ns_pirs.m, m_ns_sym1r.m, m_ns_kurch.m, m_ns_grad.m

Примітка щодо норм: у всіх солверах критерій зупинки norm(...) відповідає
евклідовій (2-й) нормі, як і в MATLAB за замовчуванням (np.linalg.norm
без параметрів теж рахує 2-норму для вектора). Норма norm(X-Xe,inf) у
головній програмі - явно нескінченна норма (np.inf), як і в оригіналі.
"""

import time
import numpy as np

# глобальний прапорець: 0 - обчислити лише F(X), 1 - обчислити (F(X), J(X))
fun = 0


# --------------------------------------------------------------------------
#  Функція правих частин і якобіан  (FdF.m, вбудована в ns_main.m)
# --------------------------------------------------------------------------
def FdF(X):
    """Обчислює f = {F(X)} (fun=0), або пару (f, J), де J - якобіан (fun=1)."""
    global fun
    N = len(X)
    m = N - 1
    f = np.zeros(N)

    f[0] = (3 + 2 * X[0]) * X[0] - 2 * X[1] - 3
    for i in range(1, m):
        f[i] = (3 + 2 * X[i]) * X[i] - X[i - 1] - 2 * X[i + 1] - 2
    f[N - 1] = (3 + 2 * X[N - 1]) * X[N - 1] - X[m - 1] - 4

    if fun == 0:
        return f

    J = np.zeros((N, N))
    J[0, 0] = 3 + 4 * X[0]
    J[0, 1] = -2
    for i in range(1, m):
        J[i, i - 1] = -1
        J[i, i] = 3 + 4 * X[i]
        J[i, i + 1] = -2
    J[N - 1, m - 1] = -1
    J[N - 1, N - 1] = 3 + 4 * X[N - 1]

    return f, J


def Dd(U, V):
    """
    Обчислює матрицю R узагальнених поділених різниць(I) 1-го порядку
    від функції F(X): R * (U - V) = F(U) - F(V).
    """
    N = len(V)
    m = N - 1
    R = np.zeros((N, N))

    z3 = 2 * U[1] + 3
    fl = (3 + 2 * U[0]) * U[0] - z3
    z2 = (3 + 2 * V[0]) * V[0]
    fr = z2 - z3
    R[0, 0] = (fl - fr) / (U[0] - V[0])

    fl = fr
    z3 = 2 * V[1] + 3
    fr = z2 - z3
    R[0, 1] = (fl - fr) / (U[1] - V[1])

    for i in range(1, m):
        im, ip = i - 1, i + 1
        z1 = 2 * U[im]
        z2 = (3 + 2 * U[i]) * U[i]
        z3 = 2 * (U[ip] + 1)
        fl = -z1 + z2 - z3
        z1 = 2 * V[im]
        fr = -z1 + z2 - z3
        R[i, im] = (fl - fr) / (U[im] - V[im])

        fl = fr
        z2 = (3 + 2 * V[i]) * V[i]
        fr = -z1 + z2 - z3
        R[i, i] = (fl - fr) / (U[i] - V[i])

        fl = fr
        z3 = 2 * (V[ip] + 1)
        fr = -z1 + z2 - z3
        R[i, ip] = (fl - fr) / (U[ip] - V[ip])

    z2 = (3 + 2 * U[N - 1]) * U[N - 1]
    fl = z2 - 4 - U[m - 1]
    z3 = 4 + V[m - 1]
    fr = z2 - z3
    R[N - 1, m - 1] = (fl - fr) / (U[m - 1] - V[m - 1])

    fl = fr
    z2 = (3 + 2 * V[N - 1]) * V[N - 1]
    fr = z2 - z3
    R[N - 1, N - 1] = (fl - fr) / (U[N - 1] - V[N - 1])

    return R


# --------------------------------------------------------------------------
#  m_ns_njut.m - метод Ньютона
# --------------------------------------------------------------------------
def m_ns_njut(f_J, X, e, mki):
    """
    Знаходження розв'язку нелiнiйної системи рiвнянь F(X) = 0
    iтерацiйним м.Ньютона:  J(X(k))*(X(k+1)-X(k)) = -F(X(k)).
    Умова закiнчення iтерацiй: norm(X(k+1)-X(k)) < e.
    """
    global fun
    fun = 1
    if e < 0 or e > 1:
        e = 1e-6
    for k in range(1, mki + 1):
        O, J = f_J(X)
        O = np.linalg.solve(J, -O)
        X = X + O
        if np.linalg.norm(O) < e:
            return X, k, 0
    return X, mki, 1


# --------------------------------------------------------------------------
#  m_ns_njutg.m - метод Ньютона з одновимірним пошуком
# --------------------------------------------------------------------------
def m_ns_njutg(f_J, X, e, mki):
    """
    Метод Ньютона + одновимірний пошук (демпфування довжини кроку).
    Умова закiнчення iтерацiй: norm(F(X(k+1))) < e.
    """
    global fun
    if e < 0 or e > 1:
        e = 1e-6
    w = None
    for k in range(1, mki + 1):
        fun = 1
        Xc = X
        F, J = f_J(X)
        P = np.linalg.solve(J, -F)
        fun = 0
        v = np.linalg.norm(F)
        for i in range(1, mki + 1):
            X = Xc + P
            w = np.linalg.norm(f_J(X))
            if w <= v:
                break
            P = 0.5 * P
        if w < e:
            return X, k, 0
    return X, mki, 1


# --------------------------------------------------------------------------
#  m_ns_broy.m - метод Бройдена (перший)
# --------------------------------------------------------------------------
def m_ns_broy(f_J, X, e, mki):
    """
    Iтерацiйний перший метод Бройдена, де X(0) - вiдоме.
    Умова закiнчення iтерацiй: norm(X(k+1)-X(k)) < e.
    """
    global fun
    fun = 1
    if e < 0 or e > 1:
        e = 1e-6
    J = None
    Dx = None
    for k in range(1, mki + 1):
        if k == 1:
            O, J = f_J(X)
            fun = 0
        else:
            O = f_J(X)
            J = J + np.outer(O, Dx) / np.dot(Dx, Dx)
        Dx = np.linalg.solve(J, -O)
        X = X + Dx
        if np.linalg.norm(Dx) < e:
            return X, k, 0
    return X, mki, 1


# --------------------------------------------------------------------------
#  m_ns_broym.m - модифікований метод Бройдена
# --------------------------------------------------------------------------
def m_ns_broym(f_J, X, e, mki):
    """
    Iтерацiйний модифікований метод Бройдена, де X(0) - вiдоме.
    Умова закiнчення iтерацiй: norm(X(k+1)-X(k)) < e.
    """
    global fun
    fun = 1
    if e < 0 or e > 1:
        e = 1e-6
    J = None
    Dx = None
    Ok = None
    for k in range(1, mki + 1):
        if k == 1:
            O, J = f_J(X)
            fun = 0
        else:
            O = f_J(X)
            G = J.T @ (O - Ok)
            J = J + np.outer(O, G) / np.dot(G, Dx)
        Dx = np.linalg.solve(J, -O)
        X = X + Dx
        if np.linalg.norm(Dx) < e:
            return X, k, 0
        Ok = O
    return X, mki, 1


# --------------------------------------------------------------------------
#  m_ns_pirs.m - метод Пірсона
# --------------------------------------------------------------------------
def m_ns_pirs(f_J, X, e, mki):
    """
    Iтерацiйний метод Пірсона, де X(0) - вiдоме.
    Умова закiнчення iтерацiй: norm(X(k+1)-X(k)) < e.
    """
    global fun
    fun = 1
    if e < 0 or e > 1:
        e = 1e-6
    J = None
    Dx = None
    Ok = None
    for k in range(1, mki + 1):
        if k == 1:
            O, J = f_J(X)
            fun = 0
        else:
            O = f_J(X)
            G = O - Ok
            J = J + np.outer(O, G) / np.dot(G, Dx)
        Dx = np.linalg.solve(J, -O)
        X = X + Dx
        if np.linalg.norm(Dx) < e:
            return X, k, 0
        Ok = O
    return X, mki, 1


# --------------------------------------------------------------------------
#  m_ns_sym1r.m - симетричний метод 1-го рангу
# --------------------------------------------------------------------------
def m_ns_sym1r(f_J, X, e, mki):
    """
    Iтерацiйний симетричний метод 1-го рангу, де X(0) - вiдоме.
    Умова закiнчення iтерацiй: norm(X(k+1)-X(k)) < e.
    """
    global fun
    fun = 1
    if e < 0 or e > 1:
        e = 1e-6
    J = None
    Dx = None
    Ok = None
    for k in range(1, mki + 1):
        if k == 1:
            O, J = f_J(X)
            fun = 0
        else:
            O = f_J(X)
            G = O - Ok - J.T @ Dx
            J = J + np.outer(O, G) / np.dot(G, Dx)
        Dx = np.linalg.solve(J, -O)
        X = X + Dx
        if np.linalg.norm(Dx) < e:
            return X, k, 0
        Ok = O
    return X, mki, 1


# --------------------------------------------------------------------------
#  m_ns_kurch.m - метод Курчатова
# --------------------------------------------------------------------------
def m_ns_kurch(F, Divdif, X, e, mki, X1):
    """
    Iтерацiйний метод Курчатова:
       F(2*X(k)-X(k-1); X(k-1)) * (X(k+1)-X(k)) = -F(X(k)),
    де X(-1), X(0) - вiдомі (X1, X), а F(U;V) - узагальнені поділені
    різниці(I) першого порядку (функція Divdif).
    Умова закiнчення iтерацiй: norm(X(k+1)-X(k)) < e.
    """
    global fun
    fun = 0
    if e < 0 or e > 1:
        e = 1e-6
    Fkm = F(X1)
    Fk = F(X)
    for k in range(1, mki + 1):
        U = 2 * X - X1
        R = Divdif(U, X1)
        O = np.linalg.solve(R, -Fk)
        X1 = X
        X = X + O
        if np.linalg.norm(O) < e:
            return X, k, 0
        Fkm = Fk
        Fk = F(X)
    return X, mki, 1


# --------------------------------------------------------------------------
#  m_ns_grad.m - градієнтний метод
# --------------------------------------------------------------------------
def m_ns_grad(f_J, X, e, mki):
    """
    Iтерацiйний градiєнтний метод, де X(0) - вiдоме.
    Розр.формули: Б.П.Демидович, И.А.Марон "Основы вычислительной
    математики", Наука, М., 1966, стор. 485-490.
    Умова закiнчення iтерацiй: norm(X(k+1)-X(k)) < e.
    """
    global fun
    fun = 1
    if e < 0 or e > 1:
        e = 1e-6
    for k in range(1, mki + 1):
        z, J = f_J(X)
        d = J.T @ z
        v = J @ d
        m = np.dot(z, v) / np.dot(v, v)
        d = m * d
        X = X - d
        if np.linalg.norm(d) < e:
            return X, k, 0
    return X, mki, 1


# --------------------------------------------------------------------------
#  Допоміжний друк і головна програма (ns_main.m)
# --------------------------------------------------------------------------
def print_result(txt, it, t, d):
    print('%-18s %8u %8.4f %20.12e' % (txt, it, t, d))


def ns_main():
    global fun

    N = 200
    Xe = np.ones(N)
    X0 = 0.5 * Xe

    e = 1e-10
    mki = 200
    pr = False  # True - друк розв'язку і значення функції

    print('Точність ітераційного процесу =' + str(e))
    print('Обмеження кількості ітерацій =' + str(mki))
    print('Кількість невідомих =' + str(N))
    print('%-18s %-8s %-8s %-20s' %
          ("Метод  розв'язання", 'Кільк.іт', 'Час(сек)', 'Похибка'))

    def show(name, X, it, t):
        print_result(name, it, t, np.linalg.norm(X - Xe, np.inf))
        if pr:
            global fun
            print(X)
            fun = 0
            print(FdF(X))

    X = X0.copy(); t1 = time.perf_counter()
    X, it, err = m_ns_njut(FdF, X, e, mki); t2 = time.perf_counter()
    show('Ньютона', X, it, t2 - t1)

    X = X0.copy(); t1 = time.perf_counter()
    X, it, err = m_ns_njutg(FdF, X, e, mki); t2 = time.perf_counter()
    show('Ньютон+однов.пошук', X, it, t2 - t1)

    X = X0.copy(); t1 = time.perf_counter()
    X, it, err = m_ns_broy(FdF, X, e, mki); t2 = time.perf_counter()
    show('Бройдена', X, it, t2 - t1)

    X = X0.copy(); t1 = time.perf_counter()
    X, it, err = m_ns_pirs(FdF, X, e, mki); t2 = time.perf_counter()
    show('Пірсона', X, it, t2 - t1)

    X = X0.copy(); t1 = time.perf_counter()
    X, it, err = m_ns_broym(FdF, X, e, mki); t2 = time.perf_counter()
    show('модифік. Бройдена', X, it, t2 - t1)

    X = X0.copy(); t1 = time.perf_counter()
    X, it, err = m_ns_sym1r(FdF, X, e, mki); t2 = time.perf_counter()
    show('симетр.1-го рангу', X, it, t2 - t1)

    X = X0.copy(); X1 = X0 - 1e-3; t1 = time.perf_counter()
    X, it, err = m_ns_kurch(FdF, Dd, X, e, mki, X1); t2 = time.perf_counter()
    show('Курчатова', X, it, t2 - t1)

    X = X0.copy(); t1 = time.perf_counter()
    X, it, err = m_ns_grad(FdF, X, e, mki); t2 = time.perf_counter()
    show('градієнтний', X, it, t2 - t1)


if __name__ == '__main__':
    ns_main()
