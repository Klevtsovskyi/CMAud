import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Функція та її точні похідні
# ============================================================

def u(x):
    return np.sin(x)


def du(x):
    return np.cos(x)


def d2u(x):
    return -np.sin(x)


# ============================================================
# Сітка
# ============================================================

X = np.linspace(0, 3, 31)

n = len(X)
n1 = n - 1

h = X[1] - X[0]
hh = 2 * h
h2 = h ** 2

U = u(X)
U1 = du(X)
U2 = d2u(X)


# Масиви для наближених похідних
Yl = np.zeros(n)       # ліва різниця
Yr = np.zeros(n)       # права різниця
Yz = np.zeros(n)       # центральна різниця
Ylr = np.zeros(n)      # друга різницева похідна
R = np.zeros(n)


# ============================================================
# Чисельні похідні
# ============================================================

for i in range(1, n - 1):

    im = i - 1
    ip = i + 1

    # Ліва різниця
    Yl[i] = (U[i] - U[im]) / h

    # Права різниця
    Yr[i] = (U[ip] - U[i]) / h

    # Центральна різниця
    Yz[i] = (U[ip] - U[im]) / hh

    # Друга різницева похідна
    Ylr[i] = (U[im] - 2 * U[i] + U[ip]) / h2


# Крайні точки
Yl[n - 1] = (U[n - 1] - U[n - 2]) / h
Yr[0] = (U[1] - U[0]) / h


# ============================================================
# Ліва різниця
# ============================================================

I = np.arange(1, n)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X[I], U1[I], "k:", label="U'(X)")
plt.plot(X[I], Yl[I], "b", label="Ліва різниця")

plt.grid()
plt.xlabel("X")
plt.ylabel("U'")
plt.title("U'(X) і ліва різниця")
plt.legend()


R[I] = U1[I] - Yl[I]

r = np.max(np.abs(R[I]))

print(f"Ліва різниця:")
print(f"h = {h:.8f}")
print(f"Теоретична похибка ~ h/2 = {h / 2:.10f}")
print(f"Розрахована похибка = {r:.10f}")


plt.subplot(1, 2, 2)

plt.plot(X[I], R[I])

plt.grid()
plt.xlabel("X")
plt.ylabel("R(X)")
plt.title("Похибка R(X)")

plt.tight_layout()
plt.show()


# ============================================================
# Права різниця
# ============================================================

I = np.arange(0, n - 1)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X[I], U1[I], "k:", label="U'(X)")
plt.plot(X[I], Yr[I], "b", label="Права різниця")

plt.grid()
plt.xlabel("X")
plt.ylabel("U'")
plt.title("U'(X) і права різниця")
plt.legend()


R[I] = U1[I] - Yr[I]

r = np.max(np.abs(R[I]))

print("\nПрава різниця:")
print(f"h = {h:.8f}")
print(f"Теоретична похибка ~ h/2 = {h / 2:.10f}")
print(f"Розрахована похибка = {r:.10f}")


plt.subplot(1, 2, 2)

plt.plot(X[I], R[I])

plt.grid()
plt.xlabel("X")
plt.ylabel("R(X)")
plt.title("Похибка R(X)")

plt.tight_layout()
plt.show()


# ============================================================
# Центральна різниця
# ============================================================

I = np.arange(1, n - 1)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X[I], U1[I], "k:", label="U'(X)")
plt.plot(X[I], Yz[I], "b", label="Центральна різниця")

plt.grid()
plt.xlabel("X")
plt.ylabel("U'")
plt.title("U'(X) і центральна різниця")
plt.legend()


R[I] = U1[I] - Yz[I]

r = np.max(np.abs(R[I]))

print("\nЦентральна різниця:")
print(f"h = {h:.8f}")
print(f"Теоретична похибка ~ h²/6 = {h2 / 6:.10f}")
print(f"Розрахована похибка = {r:.10f}")


plt.subplot(1, 2, 2)

plt.plot(X[I], R[I])

plt.grid()
plt.xlabel("X")
plt.ylabel("R(X)")
plt.title("Похибка R(X)")

plt.tight_layout()
plt.show()


# ============================================================
# Друга різниця
# ============================================================

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X[I], U2[I], "k:", label="U''(X)")
plt.plot(X[I], Ylr[I], "b", label="Друга різниця")

plt.grid()
plt.xlabel("X")
plt.ylabel("U''")
plt.title("U''(X) і друга різниця")
plt.legend()


R[I] = U2[I] - Ylr[I]

r = np.max(np.abs(R[I]))

print("\nДруга різниця:")
print(f"h = {h:.8f}")
print(f"Теоретична похибка ~ h²/12 = {h2 / 12:.10f}")
print(f"Розрахована похибка = {r:.10f}")


plt.subplot(1, 2, 2)

plt.plot(X[I], R[I])

plt.grid()
plt.xlabel("X")
plt.ylabel("R(X)")
plt.title("Похибка R(X)")

plt.tight_layout()
plt.show()


# ============================================================
# Поліноміальна апроксимація
# ============================================================

P = np.polyfit(X, U, n1)

# Перша похідна полінома
dP = np.polyder(P)

# Значення першої похідної
dY = np.polyval(dP, X)

# Друга похідна
d2P = np.polyder(dP)

# Значення другої похідної
d2Y = np.polyval(d2P, X)


# ============================================================
# Перша похідна полінома
# ============================================================

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X, U1, "k:", label="U'(X)")
plt.plot(X, dY, "b", label="dP")

plt.grid()
plt.xlabel("X")
plt.ylabel("U'")
plt.title("U'(X) і dP")
plt.legend()


R = U1 - dY

r = np.max(np.abs(R))

print("\nПоліном — перша похідна:")
print(f"Розрахована похибка = {r:.10f}")


plt.subplot(1, 2, 2)

plt.plot(X, R)

plt.grid()
plt.xlabel("X")
plt.ylabel("R(X)")
plt.title("Похибка R(X)")

plt.tight_layout()
plt.show()


# ============================================================
# Друга похідна полінома
# ============================================================

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(X, U2, "k:", label="U''(X)")
plt.plot(X, d2Y, "b", label="d2P")

plt.grid()
plt.xlabel("X")
plt.ylabel("U''")
plt.title("U''(X) і d2P")
plt.legend()


R = U2 - d2Y

r = np.max(np.abs(R))

print("\nПоліном — друга похідна:")
print(f"Розрахована похибка = {r:.10f}")


plt.subplot(1, 2, 2)

plt.plot(X, R)

plt.grid()
plt.xlabel("X")
plt.ylabel("R(X)")
plt.title("Похибка R(X)")

plt.tight_layout()
plt.show()