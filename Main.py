#finie volume code for heat equation 2 D uniform mesh by Lucka Barbeau & Matthew Coffey

import numpy as np
import get_grid
import get_heat_matrix 
import get_heat_rhs 

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt



Xmin=0
Xmax=1
Ymin=0
Ymax=1
boundary=[Xmin,Xmax,Ymin,Ymax]

nb_dof_per_dim=10

mesh=get_grid.mesh(Xmin,Xmax,Ymin,Ymax,nb_dof_per_dim)
K=np.ones((mesh.points.shape[0],1))

A=get_heat_matrix.heat_matrix(mesh,K,boundary)
B=get_heat_rhs.heat_rhs(mesh,K,boundary)
U=np.linalg.solve(A,B)



x=mesh.points[:,0]
y=mesh.points[:,1]
z=U[:,0]

fig = plt.figure()
ax = fig.gca(projection='3d')

surf=ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True,cmap=plt.cm.plasma)

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
