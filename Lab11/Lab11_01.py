"""Lab11_01: явна схема (fdm_p10), модельна задача."""
import time
import numpy as np
from heat import fdm_p10
from _lab12_common import a, b, p, f, u0, Ga, Gb, v, report

n = 101
X = np.linspace(a, b, n)
tau, M = 4e-5, 78000
J = np.arange(0, M + 1, 1200)

t1 = time.perf_counter()
U = fdm_p10(p, f, u0, Ga, Gb, v, X, J, tau)
print(f'Тривалість розрахунку ={time.perf_counter() - t1:13.3f} сек')
report(X, J, tau, U)
