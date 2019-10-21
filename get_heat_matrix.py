import numpy as np

def heat_matrix(mesh,K,boundary):
    A=np.zeros((mesh.points.shape[0],mesh.points.shape[0]))
    for i in range(mesh.points.shape[0]):
        if mesh.points[i,0]==boundary[0] or mesh.points[i,0]==boundary[1] or mesh.points[i,1]==boundary[2] or mesh.points[i,1]==boundary[3]:
            A[i,i]=1
        else:
            Aw=K[int(mesh.neighbors[i,1])]*mesh.dy/mesh.dx/(mesh.dy*mesh.dx)
            Ae=K[int(mesh.neighbors[i,0])]*mesh.dy/mesh.dx/(mesh.dy*mesh.dx)
            An=K[int(mesh.neighbors[i,3])]*mesh.dx/mesh.dy/(mesh.dy*mesh.dx)
            As=K[int(mesh.neighbors[i,2])]*mesh.dx/mesh.dy/(mesh.dy*mesh.dx)
            A[i,i]=Aw+Ae+An+As
            A[i,int(mesh.neighbors[i,0])]=-Ae
            A[i,int(mesh.neighbors[i,1])]=-Aw
            A[i,int(mesh.neighbors[i,2])]=-As
            A[i,int(mesh.neighbors[i,3])]=-An
    return A
                
                
                    
        
        