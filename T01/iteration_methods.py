import numpy as np


def jacobi(A, f, eps=1e-6):
    # print(A)
    assert np.array_equal(A, A.T)
    # assert np.greater(A, 0).all()
    # for i in range(len(A)):
    #     assert 2 * abs(A[i, i]) > np.sum(abs(A[i]))

    x = np.zeros(len(A))
    A1 = np.tril(A, -1)
    A2 = np.triu(A, 1)
    # D = np.diag(np.diag(A))
    D1 = np.diag(1 / np.diag(A))

    while True:
        x1 = -D1 @ (A1 + A2) @ x + D1 @ f
        if np.linalg.norm(x1 - x) < eps:
            break

        x = x1
    return x
