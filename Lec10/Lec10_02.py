"""Lec10_02: те саме рівняння (неявна схема) + порівняння з явною (з fileYJ.npy) і точним."""
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
tau, M = 0.01, 500
tau2 = 0.5 * tau

C = c * tau / h
B = 1 / (1 + C)
A = C * B
C = tau * B

t = 0.0
y = u0(x)
U0 = y.copy()
U1 = None
for j in range(1, M + 1):
    t += tau
    y[0] = g1(t)
    for k in range(1, N):               # неявна схема: використовує вже нове y[k-1]
        y[k] = A * y[k - 1] + B * y[k] + C * f(x[k] - h2, t - tau2)
    if 2 * j == M:
        U1 = y.copy()

t = M * tau
tm = t / 2
T1, T = f'{tm:g}', f'{t:g}'
plt.plot(x, U0, 'b.-', x, U1, 'k.--', x, y, 'r*-')
plt.xlabel('x'); plt.ylabel('u')
plt.title("Наближений розв'язок (неявна РС)")
plt.legend(['u(x,0)', f'u(x,{T1})', f'u(x,{T})'], loc='best')
plt.grid(True)
plt.show()

def u(x, t):
    return (t + 1) * x * (3 - x)        # точний розв'язок


U = u(x, t)
Yimp = y                                 # неявна схема
try:
    Yexp = np.load('fileYJ.npy')         # явна схема (запустіть спочатку le11_1.py)
except FileNotFoundError:
    raise SystemExit('Не знайдено fileYJ.npy - спочатку запустіть le11_1.py')

plt.figure()
plt.plot(x, Yexp, 'b.-', x, U, 'k.--')
plt.xlabel('x'); plt.ylabel('u')
plt.title(f"Набл.роз-к(явна РС) і точний розв'язок для t={T}")
plt.legend(['явна', 'модель'], loc='best'); plt.grid(True)
plt.show()

plt.figure()
plt.plot(x, np.abs(Yexp - U))
plt.xlabel('x'); plt.ylabel('u')
plt.title(f'Похибка |явна РС-модель| при t={T}')
plt.legend(['|явна-модель|'], loc='best'); plt.grid(True)
plt.show()
print(f'||явна-модель||_inf   = {np.max(np.abs(Yexp - U)):.4e}')

plt.figure()
plt.plot(x, Yimp, 'b.-', x, U, 'k.--')
plt.xlabel('x'); plt.ylabel('u')
plt.title(f"Набл.роз-к(неявна РС) і точний розв'язок для t={T}")
plt.legend(['неявна', 'модель'], loc='best'); plt.grid(True)
plt.show()

plt.figure()
plt.plot(x, np.abs(Yimp - U))
plt.xlabel('x'); plt.ylabel('u')
plt.title(f'Похибка |неявна РС-модель| при t={T}')
plt.legend(['|неявна-модель|'], loc='best'); plt.grid(True)
plt.show()
print(f'||неявна-модель||_inf = {np.max(np.abs(Yimp - U)):.4e}')
