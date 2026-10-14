import numpy as np
import matplotlib.pyplot as plt
from odesolvers import euler_ev, euler_nev, symm_nev


def F(x, U):
    return np.array([np.exp(-x) - U[0]])


Y0 = [1]                      # початкові умови
X = np.linspace(0, 5, 501)    # рівномірна сітка з h=0.01
Y1 = euler_ev(F, X, Y0)[0]            # явний Ейлера
epsx = epsf = 1e-12
Y2 = euler_nev(F, X, Y0, epsx, epsf)[0]   # неявний Ейлера
Y3 = symm_nev(F, X, Y0, epsx, epsf)[0]    # симетрична схема

Z = np.exp(-X) * (1 + X)      # точний розв'язок
plt.plot(X, Y1, X, Y2, X, Y3, X, Z); plt.grid(True)
plt.title("u'(x)=exp(-x)-u"); plt.xlabel('x'); plt.ylabel('u')
plt.legend(['Ейлер яв', 'Ейлер неяв', 'симет неяв', 'точне'], loc='best')
plt.show()

Z1, Z2, Z3 = Z - Y1, Z - Y2, Z - Y3
fig, ax = plt.subplots(3, 1)
for a_, z_, t in zip(ax, (Z1, Z2, Z3), ('u(x)-Y1', 'u(x)-Y2', 'u(x)-Y3')):
    a_.plot(X, z_); a_.grid(True); a_.set_title(t)
fig.tight_layout()
plt.show()

h = X[1] - X[0]
d1, d2, d3 = (np.max(np.abs(z)) for z in (Z1, Z2, Z3))
print(f"{'h':>12} {'h^2':>12} {'||u-Y1||c':>12} {'||u-Y2||c':>12} {'||u-Y3||c':>12}")
print(f'{h:12.4e} {h*h:12.4e} {d1:12.4e} {d2:12.4e} {d3:12.4e}')
