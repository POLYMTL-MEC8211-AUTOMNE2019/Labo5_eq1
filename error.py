# finie volume code for heat equation 2 D uniform mesh by Lucka Barbeau & Matthew Coffey

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 13:33:42 2019

@author: lucka
"""
import numpy as np


def error_l1(S_a, S_fv):
    error = 0
    for i in range(S_a.shape[0]):
        error = np.abs(S_a[i] - S_fv[i]) / S_a.shape[0] + error

    return error


def error_l2(S_a, S_fv):
    error = 0
    for i in range(S_a.shape[0]):
        error = ((S_a[i] - S_fv[i]) ** 2 / S_a.shape[0]) ** 0.5 + error

    return error


def error_linf(S_a, S_fv):
    error = np.max(np.abs(S_a - S_fv))
    return error
