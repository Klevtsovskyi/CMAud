from numint import i_simpn, i_trapn, i_crectn


def f(x):
    return x * abs(x)


def fm(x):
    return -x * x


def fp(x):
    return x * x


a, b = -1, 2
S = 7 / 3
Nn, Nm, Np = 512, 5, 11

print(f'Точне   = {S:.9e}')
Ss = i_simpn(a, 0, fm, Nm) + i_simpn(0, b, fp, Np)
print(f'Вар. A) = {Ss:.9e}')
print('Вар. B)')
print('  n КФ трапецій      КФ центр.прямок.')
n = 1
while n < Nn:
    n *= 2
    t = i_trapn(a, b, f, n)
    z = i_crectn(a, b, f, n)
    print(f'{n:3d} {t:.9e} {z:.9e}')
