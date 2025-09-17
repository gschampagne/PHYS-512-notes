#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 14:54:50 2025

@author: Grace
"""

import numpy as np
from matplotlib import pyplot as plt
plt.ion()

def f(x):
    return 1/x

x=np.linspace(1,2,22)
y=f(x)
ans=np.log(x[-1])-np.log(x[0]) # analytic answer for log(x)

dx=x[1]-x[0] # the spacing fo our points
myans=(0.5*(y[0]+y[-1])+np.sum(y[1:-1]))*dx

print('got:', myans, 'expected:', ans, 'error:', myans-ans)
