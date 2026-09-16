import numpy as np

g = 10
k = 15
d = 8

A = np.array([
    [3 + 1e-2 * g, 0.20, 0.40, 0.45],
    [0.20, 5 - 1e-3 * k, 0.75, 0.73],
    [0.40, 0.75, 2 - 1e-4 * d, 0.30],
    [0.45, 0.73, 0.30, 4.55]
], dtype=np.float64)

b = np.array([
    220, 100, 300, 150
], dtype=np.float64)
