import numpy as np
import scipy

from T1.iteration_methods import jacobi
from helper import *

A = A / 6
b = b / 6

E = np.eye(len(A))
A0 = E - A
# print(A0)

# 1 a
B = np.linalg.inv(A0)
# print(B)

# 1 b
# TODO (ДЗ): 2a) повинна існувати матриця 𝑩 і вона має бути невід’ємною (𝑩 ≥ 0)
assert np.linalg.norm(A, 2) < 1, "Не виконується умова збіжності, що найбільше власне число менше 1"
# TODO (ДЗ): 2d) перевірити, що усі головні мінори матриці A0 мають бути додатніми.
eps = 1e-6
S = E
P = E
while np.linalg.norm(P, 2) > eps:
    P = P @ A
    S = S + P

# print(S)

# 3 a
X = B @ b
print(X)
X = S @ b
print(X)

# 3 b I 1
X = np.linalg.solve(A0, b)
print(X)

# 3 b I 2
from scipy.linalg import lu
P, L, U = lu(A0)
y = scipy.linalg.solve_triangular(L, P @ b, lower=True)
x = scipy.linalg.solve_triangular(U, y, lower=False)
print(x)

# 3 b I 3
from scipy.linalg import qr
Q, R = qr(A0)
x = scipy.linalg.solve_triangular(R, Q.T @ b)
print(x)

# 3 b I 4
from scipy.linalg import cholesky
L = cholesky(A0, lower=True)
U = cholesky(A0, lower=False)
y = scipy.linalg.solve_triangular(L, b, lower=True)
x = scipy.linalg.solve_triangular(U, y, lower=False)
print(x)

# 3 b I 5
# TODO: виконати завдання, використовуючи sympy.linsolve() – солвер для розв’язання СЛАР.

# 3 b II 1
x = jacobi(A0, b)
print(x)
