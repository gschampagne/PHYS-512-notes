#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 09:28:15 2025

@author: Grace
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def fun(x):
    return 1/np.log(x)

x = np.linspace(2,3,21)
#print(x)
dx = x[1]-x[0]
y = fun(x)

spln = CubicSpline(x, y)
xx = np.linspace(x[0], x[-1], 2001)
yy = spln(xx)

truth = fun(xx)
plt.clf()
plt.plot(x,y,'*', label='function points')
plt.plot(xx,yy, label='spline')
plt.plot(xx,truth, label='true function')
plt.legend()
plt.show()

print("error is", np.std(truth - yy)) 
# error is only like 1 pt per million
# decr. number of points greatly increased error
