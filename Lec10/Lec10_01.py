"""Lec10_01: гіперболічне рівняння 1-го порядку u_t + c u_x = f з постійним c (явна схема).
Результат y(x, M*tau) зберігається у fileYJ.npy (його читає Lec10_02.py)."""
import numpy as np
import matplotlib.pyplot as plt

c = 1.0


def f(x, t):
    return 3 + x * (1 - x) + t * (3 - 2 * x)


def u0(x):
    return x * (3 - x)


def g1(t):
    return 0.0


a, b = 0, 3
N = 301
x = np.linspace(a, b, N)
h = x[1] - x[0]
h2 = 0.5 * h
tau, M = 0.001, 5000
tau2 = 0.5 * tau

A = c * tau / h
B = 1 - A

t = 0.0
y = u0(x)
U0 = y.copy()
U1 = None
xm = x[1:] - h2                         # вузли x(k)-h/2, k=2..N
for j in range(1, M + 1):
    t += tau
    yold = y
    y = np.empty(N)
    y[0] = g1(t)
    # явна схема: усі праві частини беруть значення з попереднього шару
    y[1:] = A * yold[:-1] + B * yold[1:] + tau * f(xm, t - tau2)
    if 2 * j == M:
        U1 = y.copy()
np.save('fileYJ.npy', y)

# графіки розв'язку при t=0, M*tau/2, M*tau
t = M * tau
tm = t / 2
plt.plot(x, U0, 'b.-', x, U1, 'k.--', x, y, 'r*-')
plt.xlabel('x'); plt.ylabel('u')
plt.title("Наближений розв'язок (явна РС)")
plt.legend(['u(x,0)', f'u(x,{tm:g})', f'u(x,{t:g})'], loc='best')
plt.grid(True)
plt.show()

U = (t + 1) * x * (3 - x)               # точний розв'язок
print(f'Похибка явної РС при t={t:g}: {np.max(np.abs(y - U)):.4e}')
