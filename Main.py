#finie volume code for heat equation 2 D uniform mesh by Lucka Barbeau & Matthew Coffey

import numpy as np
import get_grid
import get_heat_matrix 
import get_heat_rhs 
from sympy import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from eval_vector import *
from error import *

Xmin=0
Xmax=25
Ymin=0
Ymax=15
boundary=[Xmin,Xmax,Ymin,Ymax]






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

T_hat = T0 + Tx * cos(ax * pi * x / L) + Ty * sin(ay * pi * y / L) + Txy * sin(axy * pi * x * y / L**2)

k0 = 1
kx = 4
ky = 3
kxy = 5
bx = 1/2
by = 1/3
bxy = 1/4

k_hat = k0 + kx * sin(bx * pi * x / L) + ky * cos(by * pi * y / L) + kxy * cos(bxy * pi * x * y / L**2)

heat_eq_a = diff(1 * diff(T_hat, x), x) + diff(1 * diff(T_hat, y), y)

heat_eq_b = diff(k_hat * diff(T_hat, x), x) + diff(k_hat * diff(T_hat, y), y)

T_manufacturee = lambdify([x, y], T_hat, 'numpy')
k_manufacturee = lambdify([x, y], k_hat, 'numpy')

f_a = lambdify([x, y], heat_eq_a, 'numpy')
f_b = lambdify([x, y], heat_eq_b, 'numpy')










nb_dof_per_dim=60

mesh=get_grid.mesh(Xmin,Xmax,Ymin,Ymax,nb_dof_per_dim)



K_a=np.ones((mesh.points.shape[0],1))
K_b=eval_vector(k_manufacturee,mesh.points[:,0],mesh.points[:,1])


T_exact=eval_vector(T_manufacturee,mesh.points[:,0],mesh.points[:,1])
    
source_a=eval_vector(f_a,mesh.points[:,0],mesh.points[:,1])   
source_b=eval_vector(f_b,mesh.points[:,0],mesh.points[:,1])        
source_cas_base=K_a    
boundary_check=np.zeros((mesh.points.shape[0],1))   
    
    
A=get_heat_matrix.heat_matrix(mesh,K_a,boundary)
B=get_heat_rhs.heat_rhs(mesh,K_a,boundary,source_a,T_exact)
U=np.linalg.solve(A,B)



x=mesh.points[:,0]
y=mesh.points[:,1]
z=T_exact[:,0]
z2=U[:,0]

print(error_linf(z,z2))

fig = plt.figure()
ax = fig.gca(projection='3d')

surf=ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True,cmap=plt.cm.plasma)

#U_bar = plt.contourf(evalX, evalY, np.transpose(sol_rho), 10, cmap=plt.cm.plasma, origin='lower')
cbar = plt.colorbar(surf)
cbar.set_label('temp')


plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')

surf=ax.plot_trisurf(x, y, z2, linewidth=0.2, antialiased=True,cmap=plt.cm.plasma)

#U_bar = plt.contourf(evalX, evalY, np.transpose(sol_rho), 10, cmap=plt.cm.plasma, origin='lower')
cbar = plt.colorbar(surf)
cbar.set_label('temp')


plt.show()


#mesh_graph=np.meshgrid
#fig_sec = plt.figure()
#sol_sec_plot = plt.contourf(evalX, evalY, np.transpose(sol_sec), 10, cmap=plt.cm.plasma, origin='lower')
#sol_sec_plot_2 = plt.contour(sol_sec_plot, levels=sol_sec_plot.levels[::1], colors='k', origin='lower')
#cbar = plt.colorbar(sol_sec_plot)
#cbar.set_label('sec')
#cbar.add_lines(sol_sec_plot_2)
#plt.xlabel("x/l")
#plt.ylabel("y/l")
#plt.title("analytics solution sec")
#plt.show()
