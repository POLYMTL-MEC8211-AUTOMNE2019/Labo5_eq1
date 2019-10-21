# function that return a uniform grid
import numpy as np

class mesh:
    def __init__(self,Xmin,Xmax,Ymin,Ymax,nb_per_dim):
        #define points
        x = np.linspace(Xmin, Xmax, nb_per_dim)
        y = np.linspace(Ymin, Ymax, nb_per_dim)
        dx=(Xmax-Xmin)/(nb_per_dim-1)
        dy=(Ymax-Ymin)/(nb_per_dim-1)
        self.dx=dx
        self.dy=dy
        self.points=np.zeros((nb_per_dim**2,2))
        for i in range(nb_per_dim):
            for j in range(nb_per_dim):
                self.points[i*nb_per_dim+j,0]=x[i]
                self.points[i*nb_per_dim+j,1]=y[j]
                
        #define neighbors
        self.neighbors=np.zeros((nb_per_dim**2,4))
        for i in range(nb_per_dim**2):
            if self.points[i,0]!=Xmin and self.points[i,0]!=Xmax  and self.points[i,1]!=Ymin and self.points[i,1]!=Ymax:
                self.neighbors[i,0]=i+nb_per_dim
                self.neighbors[i,1]=i-nb_per_dim
                self.neighbors[i,2]=i+1
                self.neighbors[i,3]=i-1
                    

