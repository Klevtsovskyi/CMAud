import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Числове диференціювання
# ============================================================

n = 101
x = np.linspace(0, 2 * np.pi, n)

n1 = n - 1
h = x[1] - x[0]
hh = h ** 2
h2 = 2 * h

y = np.sin(x)
dy = np.cos(x)


# 1-а похідна
Dy = np.diff(y)

# Для точок x(1:n1) — права різниця
dydx = Dy / h
I = np.arange(n1)

# Центральна різниця
dydxc = (y[2:n] - y[0:n - 2]) / h2


plt.figure(figsize=(8, 8))

plt.subplot(3, 1, 1)
plt.plot(x[I], dydx - dy[I])
plt.grid()
plt.title("Похибка правої різниці")

plt.subplot(3, 1, 2)

I = np.arange(1, n1)

plt.plot(x[I], dydxc - dy[I])
plt.grid()
plt.title("Похибка центральної різниці")


# 2-а похідна
d2ydx2 = np.zeros(n - 2)

for k in range(1, n - 1):
    k1 = k - 1
    d2ydx2[k1] = (
        y[k1] - 2 * y[k] + y[k + 1]
    ) / hh


plt.subplot(3, 1, 3)

plt.plot(x[I], d2ydx2 + y[I])
plt.grid()
plt.title("Похибка 2-ї різниці")

plt.tight_layout()
plt.show()


# ============================================================
# Градієнт функції
# ============================================================

def f(x, y):
    return 16 - x ** 4 - y ** 4


x = np.linspace(-2, 2, 41)
hx = x[1] - x[0]

y = np.linspace(-2, 2, 21)
hy = y[1] - y[0]

X, Y = np.meshgrid(x, y)

Z = f(X, Y)

plt.figure()

plt.contour(Z, 10)
plt.quiver(X, Y)

plt.grid()
plt.title("f(x,y) = 16 - x^4 - y^4")

plt.show()


# ============================================================
# Поверхня
# ============================================================

fig = plt.figure()

ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, Z)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

ax.set_title("f(x,y) = 16 - x^4 - y^4")

plt.show()


# ============================================================
# Лапласіан
# ============================================================

# MATLAB:
# L = 4 .* del2(Z, hx, hy)

# Для рівномірної сітки використаємо другі похідні
d2Z_dx2 = np.gradient(
    np.gradient(Z, hx, axis=1),
    hx,
    axis=1
)

d2Z_dy2 = np.gradient(
    np.gradient(Z, hy, axis=0),
    hy,
    axis=0
)

L = d2Z_dx2 + d2Z_dy2

# MATLAB множить результат del2 на 4
L = 4 * L


fig = plt.figure()

ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, L)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("L")

ax.set_title("Δf")

plt.show()
