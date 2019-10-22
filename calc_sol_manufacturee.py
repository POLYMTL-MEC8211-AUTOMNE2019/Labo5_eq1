from sympy import *

k, T, x, y = symbols('k T x y')

# T0, Tx, Ty, Txy, ax, ay, axy, L = symbols('T0, Tx, Ty, Txy, ax, ay, axy, L')

T0 = 400
Tx = 45
Ty = 35
Txy = 27.5
ax = 1/3
ay = 1/4
axy = 1/2
L = 5

T_hat = T0 + Tx * cos(ax * pi * x / L) + Ty * sin(ay * pi * y / L) + Txy * sin(axy * pi * x * y / L)

k0 = 1
kx = 4
ky = 3
kxy = 5
bx = 1/2
by = 1/3
bxy = 1/4

k_hat = k0 + kx * sin(bx * pi * x / L) + ky * cos(by * pi * y / L) + kxy * cos(bxy * pi * x * y / L)

heat_eq_a = diff(1 * diff(T_hat, x), x) + diff(1 * diff(T_hat, y), y)

heat_eq_b = diff(k_hat * diff(T_hat, x), x) + diff(k_hat * diff(T_hat, y), y)

T_manufacturee = lambdify([x, y], T_hat, 'numpy')
k_manufacturee = lambdify([x, y], k_hat, 'numpy')

f_a = lambdify([x, y], heat_eq_a, 'numpy')
f_b = lambdify([x, y], heat_eq_b, 'numpy')

