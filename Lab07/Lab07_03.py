import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Задача Коші: u'' + 3u' + 7u = sin(u);  u(0)=3, u'(0)=-1
def f1(u1, u2):
    return u2 + 0 * u1


def f2(u1, u2):
    return np.sin(u1) - 3 * u2 - 7 * u1


def f(t, U):
    return [f1(U[0], U[1]), f2(U[0], U[1])]


U0 = [3, -1]
sol = solve_ivp(f, [0, 10], U0, method='RK45', t_eval=np.linspace(0, 10, 600))
T, u, du = sol.t, sol.y[0], sol.y[1]

# --- графіки
plt.figure(figsize=(10, 8))
plt.subplot(2, 2, 1)
plt.plot(T, u, 'b', linewidth=2, label='u(t)')
plt.plot(T, du, 'r', label="u'(t)")
plt.grid(True); plt.legend(loc='best')
plt.title("u''+3u'+7u=sin(u); u(0)=3, u'(0)=-1")
plt.xlabel('t'); plt.ylabel("u(t),u'(t)")

plt.subplot(2, 2, 2)
plt.plot(u, du, 'k', linewidth=2)
plt.axis('equal'); plt.grid(True)
plt.title('Фазова траєкторія'); plt.xlabel('u'); plt.ylabel("u'")

x = np.linspace(u.min(), u.max(), 11)
y = np.linspace(du.min(), du.max(), 11)
X, Y = np.meshgrid(x, y)
F1, F2 = f1(X, Y), f2(X, Y)

plt.subplot(2, 2, 3)
plt.quiver(X, Y, F1, F2)
plt.axis('equal'); plt.title('Поле напрямків')

plt.subplot(2, 2, 4)
plt.quiver(X, Y, F1, F2)
plt.plot(u, du, 'k', linewidth=2)
plt.axis('equal'); plt.grid(True)
plt.title('Фазова траєкторія+Поле напрямків')
plt.xlabel('u'); plt.ylabel("u'")
plt.tight_layout()
plt.show()
