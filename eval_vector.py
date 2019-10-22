# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 14:22:39 2019

@author: lucka
"""
import numpy as np

def eval_vector(f,x,y):
    z=np.zeros((x.shape[0],1))
    for i in range(x.shape[0]):
        z[i]=f(x[i],y[i])
    return z