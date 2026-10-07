import numpy as np
import matplotlib.pyplot as plt
from inteq import fr2_mzjd

# a) модельна задача: K(x,s) = x*s^2 + x^2*s
def k_a(x):
    t = x * x
    return np.array([x, t]), np.array([t, x])   # A_k(x), B_k(s)


def f_a(x):
    return (11 / 12 - x / 9) * x


lamb, a, b, n = 1 / 3, 0, 1, 21
X, Y, n = fr2_mzjd(lamb, a, b, 1, n, k_a, f_a)
plt.plot(X, Y, 'm', linewidth=2)
plt.title("Наближений розв'язок IP y")
plt.xlabel('x'); plt.ylabel('y'); plt.grid(True)
plt.show()

Z = X - Y                     # точний розв'язок u(x)=x
print(f'a) max|u-y| = {np.max(np.abs(Z)):.3e}')
plt.plot(X, Z, 'r', linewidth=2)
plt.title("Різниця розв'язків u-y")
plt.xlabel('x'); plt.ylabel('u-y'); plt.grid(True)
plt.show()


# б) вироджене ядро K(x,s) = sum_{i=1}^m x^i s^i / i!
def make_k(m):
    def k(x):
        A = np.zeros(m)
        B = np.zeros(m)
        z = 1.0
        fac = 1.0
        for i in range(1, m + 1):
            z *= x
            fac *= i
            A[i - 1] = z
            B[i - 1] = z / fac
        return A, B
    return k


def f_b(x):
    return x


lamb, a, b = 1, 0, 2
colors = 'mbkr'
labels = []
for j, m in enumerate([4, 8, 12, 16]):
    X, Y, n = fr2_mzjd(lamb, a, b, 1, n, make_k(m), f_b)
    plt.plot(X, Y, colors[j])
    labels.append(f'm={m}')
plt.title("Наближений розв'язок IP")
plt.grid(True); plt.xlabel('x'); plt.ylabel('y')
plt.legend(labels, loc='best')
plt.show()
