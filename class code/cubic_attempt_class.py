#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 09:01:37 2025

@author: Grace
"""

import numpy as np
import matplotlib.pyplot as plt

def fun(x):
    return 1/np.log(x)

x = np.linspace(2,3,21)
#print(x)
y = fun(x)

plt.clf()
plt.plot(x,y,'*')
plt.show()

xx = np.linspace(x[1], x[-2]-1e-6, 1001)
yy = 0 * xx

for i in range(len(xx)):
    # trying to see hich values to use in fit
    ind = np.max(np.where(xx[i]>=x)[0]) # find index of my left-hand neighbour
    fitp = np.polyfit(x[ind-1:ind+3], y[ind-1:ind+3], 3) # fit polynomial to left neighbour minus 1 and right neighbour +3
    yy[i] = np.polyval(fitp,xx[i])

truth = fun(xx)
print('average error: ', np.std(yy-truth))
