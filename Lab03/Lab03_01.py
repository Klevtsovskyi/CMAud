import numpy as np
import matplotlib.pyplot as plt


def u(t):
    return np.exp(np.sin(np.pi * t)) / t**3


# Параметри
n = 11
nm = n - 1
step = 5
N = nm * step + 1

a = 1.5
b = 3.5


# Рівномірна сітка
Xr = np.linspace(a, b, n)

# Нерівномірна сітка
Xn = np.sort(
    Xr[0] + (Xr[-1] - Xr[0]) * np.random.rand(n)
)

# Значення функції у вузлах
Ur = u(Xr)
Un = u(Xn)

# Згущена сітка для побудови графіка
x = np.linspace(a, b, N)

# Індекси вузлів
I = np.arange(0, N, step)


# ==========================================
# Інтерполяція на рівномірній сітці
# ==========================================

Pr = np.polyfit(Xr, Ur, nm)
Yr = np.polyval(Pr, x)


# ==========================================
# Інтерполяція на нерівномірній сітці
# ==========================================

Pn = np.polyfit(Xn, Un, nm)
Yn = np.polyval(Pn, x)


# ==========================================
# Рівномірна сітка
# ==========================================

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(Xr, Ur, 'ko:', label='Вузли')
plt.plot(x, Yr, 'b', label='Інтерполяційний поліном')

plt.grid()
plt.xlabel('Xr, x')
plt.ylabel('Ur, Yr')
plt.title('Рівномірна сітка')
plt.legend()


# Похибка
Rr = Yr[I] - Ur

plt.subplot(1, 2, 2)

plt.plot(Xr, Rr)

plt.grid()
plt.xlabel('Xr')
plt.ylabel('Rr')
plt.title('Похибка Rr = Yr - Ur')

plt.tight_layout()
plt.show()


# ==========================================
# Нерівномірна сітка
# ==========================================

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)

plt.plot(Xn, Un, 'ko:', label='Вузли')
plt.plot(x, Yn, 'b', label='Інтерполяційний поліном')

plt.grid()
plt.xlabel('Xn, x')
plt.ylabel('Un, Yn')
plt.title('Нерівномірна сітка')
plt.legend()


# Похибка
Rn = Yn[I] - Un

plt.subplot(1, 2, 2)

plt.plot(Xn, Rn)

plt.grid()
plt.xlabel('Xn')
plt.ylabel('Rn')
plt.title('Похибка Rn = Yn - Un')

plt.tight_layout()
plt.show()


# ==========================================
# Порівняння похибок
# ==========================================

plt.figure()

plt.plot(Xr, Rr, 'k', label='Рівномірна')
plt.plot(Xn, Rn, 'b', label='Нерівномірна')

plt.grid()
plt.xlabel('Xr, Xn')
plt.ylabel('Rr, Rn')
plt.title('Похибки для різних сіток')
plt.legend()

plt.show()
