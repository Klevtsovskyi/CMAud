import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline


# ============================================================
# Визначення функції та її похідних
# ============================================================

def u(x):
    return np.sin(x)


def du(x):
    return np.cos(x)


def d2u(x):
    return -np.sin(x)


def d3u(x):
    return -np.cos(x)


# ============================================================
# Задаємо сплайнову сітку X
# ============================================================

n = 16

X = np.linspace(0, 6, n)


# Значення функції та її похідних
U = u(X)
U1 = du(X)
U2 = d2u(X)
U3 = d3u(X)


# ============================================================
# Кубічний сплайн
# ============================================================

S31 = CubicSpline(X, U)


# ============================================================
# Інтерполяційна сітка
# ============================================================

N = 10 * n + 1

XX = np.linspace(0, 6, N)


# ============================================================
# u(X) та S31(XX)
# ============================================================

Y = S31(XX)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X, U, "k:", label="u(X)")
plt.plot(XX, Y, "b", label="S31(XX)")

plt.grid()
plt.xlabel("x")
plt.ylabel("u")
plt.title("u(X) і S31(XX)")
plt.legend()


# Похибка
R = Y - u(XX)

r = np.max(np.abs(R))

print(f"||S31(x;u) - u||∞ = {r:.12f}")


plt.subplot(1, 2, 2)

plt.plot(XX, R)

plt.grid()
plt.xlabel("x")
plt.ylabel("R(x)")
plt.title("Похибка R(X) для u")

plt.tight_layout()
plt.show()


# ============================================================
# Перша похідна
# ============================================================

D = S31.derivative(1)

Y = D(XX)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X, U1, "k:", label="u'(X)")
plt.plot(XX, Y, "b", label="S31'(XX)")

plt.grid()
plt.xlabel("x")
plt.ylabel("u'")
plt.title("u'(X) і S31'(XX)")
plt.legend()


R = Y - du(XX)

r = np.max(np.abs(R))

print(f"||S31'(x;u) - u'||∞ = {r:.12f}")


plt.subplot(1, 2, 2)

plt.plot(XX, R)

plt.grid()
plt.xlabel("x")
plt.ylabel("R(x)")
plt.title("Похибка R(X) для u'")

plt.tight_layout()
plt.show()


# ============================================================
# Друга похідна
# ============================================================

D = S31.derivative(2)

Y = D(XX)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X, U2, "k:", label="u''(X)")
plt.plot(XX, Y, "b", label="S31''(XX)")

plt.grid()
plt.xlabel("x")
plt.ylabel("u''")
plt.title("u''(X) і S31''(XX)")
plt.legend()


R = Y - d2u(XX)

r = np.max(np.abs(R))

print(f"||S31''(x;u) - u''||∞ = {r:.12f}")


plt.subplot(1, 2, 2)

plt.plot(XX, R)

plt.grid()
plt.xlabel("x")
plt.ylabel("R(x)")
plt.title("Похибка R(X) для u''")

plt.tight_layout()
plt.show()


# ============================================================
# Третя похідна
# ============================================================

D = S31.derivative(3)

Y = D(XX)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X, U3, "k:", label="u'''(X)")
plt.plot(XX, Y, "b", label="S31'''(XX)")

plt.grid()
plt.xlabel("x")
plt.ylabel("u'''")
plt.title("u'''(X) і S31'''(XX)")
plt.legend()


R = Y - d3u(XX)

r = np.max(np.abs(R))

print(f"||S31'''(x;u) - u'''||∞ = {r:.12f}")


plt.subplot(1, 2, 2)

plt.plot(XX, R)

plt.grid()
plt.xlabel("x")
plt.ylabel("R(x)")
plt.title("Похибка R(X) для u'''")

plt.tight_layout()
plt.show()
