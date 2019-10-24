# finie volume code for heat equation 2 D uniform mesh by Lucka Barbeau & Matthew Coffey

import numpy as np
import get_grid
import get_heat_matrix
import get_heat_rhs
from sympy import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from eval_vector import *
from error import *

Xmin = 0
Xmax = 25
Ymin = 0
Ymax = 15
boundary = [Xmin, Xmax, Ymin, Ymax]

k, T, x, y = symbols('k T x y')

# T0, Tx, Ty, Txy, ax, ay, axy, L = symbols('T0, Tx, Ty, Txy, ax, ay, axy, L')

T0 = 400
Tx = 45
Ty = 35
Txy = 27.5
ax = 1 / 3
ay = 1 / 4
axy = 1 / 2
L = 5

T_hat = T0 + Tx * cos(ax * pi * x / L) + Ty * sin(ay * pi * y / L) + Txy * sin(axy * pi * x * y / L ** 2)

k0 = 1
kx = 4
ky = 3
kxy = 5
bx = 1 / 2
by = 1 / 3
bxy = 1 / 4

k_hat = k0 + kx * sin(bx * pi * x / L) + ky * cos(by * pi * y / L) + kxy * cos(bxy * pi * x * y / L ** 2)

heat_eq_a = diff(1 * diff(T_hat, x), x) + diff(1 * diff(T_hat, y), y)

heat_eq_b = diff(k_hat * diff(T_hat, x), x) + diff(k_hat * diff(T_hat, y), y)

T_manufacturee = lambdify([x, y], T_hat, 'numpy')
k_manufacturee = lambdify([x, y], k_hat, 'numpy')

f_a = lambdify([x, y], heat_eq_a, 'numpy')
f_b = lambdify([x, y], heat_eq_b, 'numpy')

h1 = 17
mesh1 = get_grid.mesh(Xmin, Xmax, Ymin, Ymax, h1)
K_a1 = np.ones((mesh1.points.shape[0], 1))
K_b1 = eval_vector(k_manufacturee, mesh1.points[:, 0], mesh1.points[:, 1])
T_exact1 = eval_vector(T_manufacturee, mesh1.points[:, 0], mesh1.points[:, 1])
source_a1 = eval_vector(f_a, mesh1.points[:, 0], mesh1.points[:, 1])
source_b1 = eval_vector(f_b, mesh1.points[:, 0], mesh1.points[:, 1])
A1 = get_heat_matrix.heat_matrix(mesh1, K_a1, boundary)
B1 = get_heat_rhs.heat_rhs(mesh1, K_a1, boundary, source_a1, T_exact1)
U1 = np.linalg.solve(A1, B1)
x1 = mesh1.points[:, 0]
y1 = mesh1.points[:, 1]
sol_exact1 = T_exact1[:, 0]
sol_calc1 = U1[:, 0]

# print(error_linf(z, z2))
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# surf = ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True, cmap=plt.cm.plasma)
#
# # U_bar = plt.contourf(evalX, evalY, np.transpose(sol_rho), 10, cmap=plt.cm.plasma, origin='lower')
# cbar = plt.colorbar(surf)
# cbar.set_label('temp')
#
# plt.show()
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# surf = ax.plot_trisurf(x, y, z2, linewidth=0.2, antialiased=True, cmap=plt.cm.plasma)
#
# # U_bar = plt.contourf(evalX, evalY, np.transpose(sol_rho), 10, cmap=plt.cm.plasma, origin='lower')
# cbar = plt.colorbar(surf)
# cbar.set_label('temp')
#
# plt.show()

# mesh_graph=np.meshgrid
# fig_sec = plt.figure()
# sol_sec_plot = plt.contourf(evalX, evalY, np.transpose(sol_sec), 10, cmap=plt.cm.plasma, origin='lower')
# sol_sec_plot_2 = plt.contour(sol_sec_plot, levels=sol_sec_plot.levels[::1], colors='k', origin='lower')
# cbar = plt.colorbar(sol_sec_plot)
# cbar.set_label('sec')
# cbar.add_lines(sol_sec_plot_2)
# plt.xlabel("x/l")
# plt.ylabel("y/l")
# plt.title("analytics solution sec")
# plt.show()
