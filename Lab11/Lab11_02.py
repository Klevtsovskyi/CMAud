"""Lab11_02: неявна схема з вагою sigm=0.5 (fdm_p1s), модельна задача."""
import time
import numpy as np
from heat import fdm_p1s
from _lab12_common import a, b, p, f, u0, Ga, Gb, v, report

n = 101
X = np.linspace(a, b, n)
tau, M = 5e-3, 630
J = np.arange(0, M + 1, 10)

t1 = time.perf_counter()
U = fdm_p1s(p, f, u0, Ga, Gb, v, X, J, tau, 0.5)
print(f'Тривалість розрахунку ={time.perf_counter() - t1:13.3f} сек')
report(X, J, tau, U)
