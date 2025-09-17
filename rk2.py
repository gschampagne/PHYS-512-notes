#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 13:24:11 2025

@author: Grace
"""

import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    # dy/dx = f(x,y)
    # if f=-y, then dy/dx = -y, then y=exp(-x)
    return -y

def step(f,x,y,h):
    k1 = h * f(x,y)
    #k2 = h * f(x+h, y+k1)
    #return y+0.5*(k1+k2)
    k2 = h* f(x+h/2, y+k1/2) # other version using half a step
    return y+k2

n = 40
h = 1/n
x = np.linspace(0,1,n+1)
y = 0*x
y[0] = 1
for i in range(len(x)-1):
    y[i+1] = step(f,x[i],y[i],h)
    # new y is old y plus a step
plt.figure(0)
plt.clf()
plt.plot(x,y)
plt.plot(x,np.exp(-x))
plt.show()

plt.figure(1)
plt.plot(x,y-np.exp(-x))
plt.show()