#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 09:51:15 2025

@author: Grace
"""


import numpy as np
import matplotlib.pyplot as plt
plt.ion()

np.random.seed(5) # makes sure have same randomization each time

x = np.linspace(-3,3,3001)
y_true = (x**4 + 3*x**3 - 3*x**2 + 2.5*x -1) / 20

plt.clf()
plt.plot(x,y_true)
plt.title('Polynomial')
plt.show()

y = y_true + np.random.randn(len(x))
# where if we have all sigma = 1 then N becomes the identity matrix and goes away
# if noise is teh same everywhere can just ignore N
plt.plot(x,y,'.')

order = 5 # say is a 4th order polynomial
A = np.zeros([len(y), 5])
# want A*m to give me polynomial in x
# so A[:,0] = 1
# A[:,1] = x
# A[:,2] = x**2

# A * (x**0 coeff, x**1 coeff, x**2 coeff, ...)
for i in range(order):
    A[:,i] = x**i

ptrue = np.asarray([-1,2.5,-3,3,1]) / 20
yy = A@ptrue
plt.plot(x,yy)

lhs = A.T@A
rhs = A.T@y
u,s,v, = np.linalg.svd(lhs,0)
sinv = 1/s
mask = s<1e-12*s.max() # find smallish singular values
sinv[mask] = 0 # zap them out
lhs_inv = v.T@(np.diag(sinv))@u.T

#fitp = np.linalg.inv(lhs)@rhs
fitp = lhs_inv@rhs
print('ptrue: ', ptrue)
print('fitp:, ', fitp)

y_fit = A@fitp # minimize chi^2
plt.plot(x,y_fit)
print('chi^2 is:', np.sum((y_fit-y)**2))


plt.title('Polynomial Linear Least Squares + random')
plt.show()