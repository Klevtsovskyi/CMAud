"""Крайові задачі для ЗДР 2-го порядку методом сіток + методи прогонки.
(Переклад з MATLAB: m_progm.m, m_progn.m, m_grid2.m, m_grid4.m)

Рівняння  u''(x) - p(x) u(x) = -f(x),  x in (a,b),  p >= 0, з крайовими умовами
    D[0]*u'(a) + D[1]*u(a) = -D[2],
    D[3]*u'(b) + D[4]*u(b) =  D[5].
"""
import numpy as np


def _check_args(D, a, b):
    if len(D) < 6:
        raise ValueError('Недостатньо елементів у масиві D (потрібно 6)')
    if abs(D[0]) + abs(D[1]) == 0 or abs(D[3]) + abs(D[4]) == 0:
        raise ValueError('Некоректні коефіцієнти крайових умов D')
    if a >= b:
        raise ValueError('a >= b')


# ------------------------------------------------------------ прогонки
def m_progm(a, c, b, f):
    """Монотонна (права) прогонка для
        a[k]*y[k-1] - c[k]*y[k] + b[k]*y[k+1] = -f[k],   a[0]=b[n-1]=0,
    за умови |c[k]| >= |a[k]|+|b[k]|.
    Повертає (y, alfa, err); err=True - умови не виконані, розв'язок не знайдено."""
    a, c, b, f = (np.asarray(v, dtype=float) for v in (a, c, b, f))
    n = len(c)
    alfa = np.zeros(n)
    y = np.zeros(n)
    if np.any(np.abs(c) < np.abs(a) + np.abs(b)):
        return y, alfa, True
    beta = np.zeros(n)
    for k in range(n - 1):                       # пряма прогонка
        z = c[k] - a[k] * alfa[k]
        alfa[k + 1] = b[k] / z
        beta[k + 1] = (f[k] + a[k] * beta[k]) / z
    y[n - 1] = (f[n - 1] + a[n - 1] * beta[n - 1]) / (c[n - 1] - a[n - 1] * alfa[n - 1])
    for k in range(n - 2, -1, -1):               # зворотна прогонка
        y[k] = alfa[k + 1] * y[k + 1] + beta[k + 1]
    return y, alfa, False


def m_progn(a, c, b, f):
    """Немонотонна прогонка (метод Гауса з вибором провідного елемента по рядку).
    Та сама СЛАР, що й у m_progm. Повертає (y, alfa)."""
    a, c, b, f = (np.asarray(v, dtype=float) for v in (a, c, b, f))
    n = len(c)
    nm = n - 1                                   # кількість кроків
    alfa = np.zeros(n)
    beta = np.zeros(n)
    y = np.zeros(n)
    ka = np.zeros(n, dtype=int)                  # ka[0]=0
    te = np.zeros(n, dtype=int)
    C, A, F, Q = c[0], a[1], f[0], f[1]
    for k in range(nm):
        j = k + 1
        if abs(C) >= abs(b[k]):
            alfa[j] = b[k] / C
            beta[j] = F / C
            C = c[j] - A * alfa[j]
            F = Q + A * beta[j]
            te[j] = ka[k]
            ka[j] = j
            if k != nm - 1:
                A, Q = a[j + 1], f[j + 1]
        else:
            alfa[j] = C / b[k]
            beta[j] = -F / b[k]
            C = c[j] * alfa[j] - A
            F = Q - c[j] * beta[j]
            te[j] = j
            ka[j] = ka[k]
            if k != nm - 1:
                A = a[j + 1] * alfa[j]
                Q = f[j + 1] + a[j + 1] * beta[j]
    y[ka[n - 1]] = F / C
    for k in range(nm - 1, -1, -1):
        j = k + 1
        y[te[j]] = alfa[j] * y[ka[j]] + beta[j]
    return y, alfa


# ------------------------------------------------------------ різницеві схеми
def m_grid2(pfun, ffun, a, b, D, n):
    """Різницева схема 2-го порядку точності. Повертає (X, Y)."""
    _check_args(D, a, b)
    X = np.linspace(a, b, n)
    h = X[1] - X[0]
    hh, h2 = h * h, 2 * h
    A = np.ones(n)
    B = np.ones(n)
    C = np.ones(n)
    F = np.ones(n)
    A[0] = 0
    B[0] = 2 * D[0]
    x = D[0] * hh
    C[0] = B[0] - D[1] * h2 + x * pfun(a)
    F[0] = D[2] * h2 + x * ffun(a)
    for k in range(1, n - 1):
        x = X[k]
        C[k] = 2 + hh * pfun(x)
        F[k] = hh * ffun(x)
    A[n - 1] = 2 * D[3]
    B[n - 1] = 0
    x = D[3] * hh
    C[n - 1] = A[n - 1] + D[4] * h2 + x * pfun(b)
    F[n - 1] = D[5] * h2 + x * ffun(b)
    # монотонна прогонка, а за невиконання умов стійкості - немонотонна
    Y, _, err = m_progm(A, C, B, F)
    if err:
        Y, _ = m_progn(A, C, B, F)
    return X, Y


def m_grid4(p, ffun, dffun, a, b, D, n):
    """Різницева схема 4-го порядку точності для p=const>0.
    ffun(x) -> (f(x), f''(x));  dffun(x) -> f'(x). Повертає (X, Y)."""
    _check_args(D, a, b)
    X = np.linspace(a, b, n)
    h = X[1] - X[0]
    hh = h * h
    hh6, hh12 = hh / 6, hh / 12
    h4, h2 = 0.25 * h, 0.5 * h
    ph6, ph12 = 1 + p * hh6, 1 + p * hh12
    e1 = ph12 * hh
    hv = hh12 * hh
    x = 2 + p * e1
    A = np.ones(n)
    B = np.ones(n)
    C = x * np.ones(n)
    F = np.ones(n)
    # x = a
    A[0] = 0
    B[0] = D[0] / h
    C[0] = B[0] + h2 * p * D[0] * ph12 - D[1] * ph6
    f, f2 = ffun(a)
    f1 = dffun(a)
    F[0] = D[2] * ph6 + D[0] * (h2 * ph12 * f + hh6 * (f1 + h4 * f2))
    # внутрішні вузли
    for i in range(1, n - 1):
        f, f2 = ffun(X[i])
        F[i] = e1 * f + hv * f2
    # x = b
    A[n - 1] = D[3] / h
    B[n - 1] = 0
    C[n - 1] = A[n - 1] + h2 * p * ph12 * D[3] + ph6 * D[4]
    f, f2 = ffun(b)
    f1 = dffun(b)
    F[n - 1] = ph6 * D[5] + D[3] * (h2 * ph12 * f - hh6 * (f1 - h4 * f2))
    Y, _ = m_progn(A, C, B, F)
    return X, Y
