# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:29:05 2019

@author: lucka
"""
import numpy as np

def heat_rhs(mesh,rhs,boundary):
    B=np.zeros((mesh.points.shape[0],1))
    for i in range(mesh.points.shape[0]):
        if mesh.points[i,0]==boundary[0] or mesh.points[i,0]==boundary[1] or mesh.points[i,1]==boundary[2] or mesh.points[i,1]==boundary[3]:
            B[i]=0
        else:
            B[i]=1
    return B