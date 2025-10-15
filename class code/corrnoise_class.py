#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 22:24:56 2025

@author: Grace
"""

import numpy as np
import matplotlib.pyplot as plt
plt.ion()

ndata = 1000
N = np.zeros([ndata, ndata])

# we'll decrease that N_ij = exp(-0.5*(i-j)^2/width^2)
width = 30
for i in range(ndata):
    for j in range(ndata):
        N[i,j] = np.exp(-0.5*(i-j)**2/width**2)
        # matrix is symmetric
plt.figure(1)
plt.clf()
plt.imshow(N)
plt.show()
# since symmetric need to use eigh (runs twice as fats and has better answers)
e, v = np.linalg.eigh(N)
# some e are negative due to roundoff so set to zero
e[e<0] = 0


d_uncorr = np.sqrt(e)*np.random.randn(len(e)) # this uncorrelated noisy data
d = v@d_uncorr

plt.figure(2)
plt.clf()
plt.plot(d)
plt.show()

plt.clf()
plt.plot(np.sqrt(e), '.') 
plt.show()

plt.clf()
plt.plot(v[:,-3]) # eigen vector looks like chebyshev
plt.show()