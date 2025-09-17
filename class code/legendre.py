#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 15:53:02 2025

@author: Grace
"""

import numpy as np
import matplotlib.pyplot as plt

plt.ion()

x=np.linspace(-1, 1, 2001)

ord=5
P=np.zeros([len(x),ord+1])
P[:,0]=1.0 
P[:,1]=x 
for n in range(1,ord):
    P[:,n+1]=((2*n+1)*x*P[:,n]-n*P[:,n-1])/(n+1)

plt.clf()
plt.plot(x,P)
plt.show()
#plt.plot(x,P[:,2]) # Quadratic

overlap=P.T@P/len(x) 
# sum for each column of P against each other column of P 
# rougly zero on the diagonal
print(overlap)