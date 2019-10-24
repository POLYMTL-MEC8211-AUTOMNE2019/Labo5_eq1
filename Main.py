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
import math

Xmin = 0
Xmax = 5
Ymin = 0
Ymax = 3
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
L2= 3
#### IMPORTANT NOTE !!! !!!!! POUR AVOIR LA BONNE SURFACE DE REPONSE COMPARATIVEMENT AU LIVRE L DOIT ETRE ENLEVÉ  SI L'ON GARDE L ON DOIT AVOIR X DE 0 A 25 et X de 0 A 15 !!!!
T_hat = T0 + Tx * cos(ax * pi * x ) + Ty * sin(ay * pi * y  ) + Txy * sin(axy * pi * x * y  )

k0 = 1
kx = 4
ky = 3
kxy = 5
bx = 1 / 2
by = 1 / 3
bxy = 1 / 4

#### IMPORTANT NOTE !!! !!!!! POUR AVOIR UNE SOLUTION PHYSIQUE ON DOIT GARER L DANS L'éQUATION CAR SINON K DEVIENT NEGATIF , SOLUTION DIVERENTE DE L'EQUATION DIFFERENTIEL !!!!
k_hat = k0 + kx * sin(bx * pi * x / L) + ky * cos(by * pi * y / L2) + kxy * cos(bxy * pi * x * y / L )

heat_eq_a = diff(1 * diff(T_hat, x), x) + diff(1 * diff(T_hat, y), y)

heat_eq_b = diff(k_hat * diff(T_hat, x), x) + diff(k_hat * diff(T_hat, y), y)

T_manufacturee = lambdify([x, y], T_hat, 'numpy')
k_manufacturee = lambdify([x, y], k_hat, 'numpy')

f_a = lambdify([x, y], heat_eq_a, 'numpy')
f_b = lambdify([x, y], heat_eq_b, 'numpy')

mesh_sizes = np.array([17, 33, 65, 129])

error = np.zeros((2, len(mesh_sizes)))
i = 0

for h in mesh_sizes:
    mesh = get_grid.mesh(Xmin, Xmax, Ymin, Ymax, h)
    K_a = np.ones((mesh.points.shape[0], 1))
    K_b = eval_vector(k_manufacturee, mesh.points[:, 0], mesh.points[:, 1])
    T_exact = eval_vector(T_manufacturee, mesh.points[:, 0], mesh.points[:, 1])
    source_a = eval_vector(f_a, mesh.points[:, 0], mesh.points[:, 1])
    source_b = eval_vector(f_b, mesh.points[:, 0], mesh.points[:, 1])
    A_a = get_heat_matrix.heat_matrix(mesh, K_a, boundary)
    B_a = get_heat_rhs.heat_rhs(mesh, K_a, boundary, source_a, T_exact)
    A_b = get_heat_matrix.heat_matrix(mesh, K_b, boundary)
    B_b = get_heat_rhs.heat_rhs(mesh, K_b, boundary, source_b, T_exact)
    U_a = np.linalg.solve(A_a, B_a)
    U_b = np.linalg.solve(A_b, B_b)
    x1 = mesh.points[:, 0]
    y1 = mesh.points[:, 1]
    sol_exact = T_exact[:, 0]
    sol_calc_a = U_a[:, 0]
    sol_calc_b = U_b[:, 0]
    error[0, i] = error_l1(sol_exact, sol_calc_a)
    error[1, i] = error_l1(sol_exact, sol_calc_b)
    i += 1

fig2= plt.figure()
plt.loglog(mesh_sizes, error[0, :])
plt.loglog(mesh_sizes, error[1, :])
plt.xlabel("h")
plt.ylabel("Erreur")
plt.legend(("Erreur pour une conductivité thermique constante", "Erreur pour une conductivité thermique variable"))
plt.title('Ordre de convergence des erreurs')
plt.show()

##### CALCUL  ERREUR COMME 5.3.2 
p_const = math.log(error[0, -2]/error[0, -1])/math.log(mesh_sizes[-2]/mesh_sizes[-1])
p_var = math.log(error[1, -2]/error[1, -1])/math.log(mesh_sizes[-2]/mesh_sizes[-1])

print("L'ordre de convergence pour une conductivité thermique constante est : {:5.2f}".format(p_const))
print("L'ordre de convergence pour une conductivité thermique variable est : {:5.2f}".format(p_var))

x=mesh.points[:,0]
y=mesh.points[:,1]
z=sol_calc_b


#### REPRESENTATION DE LA SURFACE 
fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True, cmap=plt.cm.plasma)


cbar = plt.colorbar(surf)
cbar.set_label('temp')
plt.xlabel("x")
plt.ylabel("Y")
plt.title('Surface de reponse pour la temperature selon la solution manifacturer.')
plt.show()

