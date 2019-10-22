#finie volume code for heat equation 2 D uniform mesh by Lucka Barbeau & Matthew Coffey

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:29:05 2019

@author: lucka
"""
import numpy as np

def heat_rhs(mesh,rhs,boundary,termsource,boundary_value):
    B=np.zeros((mesh.points.shape[0],1))
    for i in range(mesh.points.shape[0]):
        if mesh.points[i,0]==boundary[0] or mesh.points[i,0]==boundary[1] or mesh.points[i,1]==boundary[2] or mesh.points[i,1]==boundary[3]:
            #PUT BOUNDARY CONDITION HERE !
            B[i]=boundary_value[i]
        else:
            #PUT SOURCE TERM HERE !
            B[i]=-termsource[i]
    return B