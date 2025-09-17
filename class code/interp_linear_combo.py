#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 15:13:19 2025

@author: Grace
"""

import numpy as np
from matplotlib import pyplot as plt
plt.ion()

def f(x):
    return 1/x

x=np.linspace(1,2,11)
y=f(x)
ans=np.log(x[-1])-np.log(x[0]) # analytic answer for log(x)

dx=x[1]-x[0] # the spacing fo our points
myans=(0.5*(y[0]+y[-1])+np.sum(y[1:-1]))*dx

print('got:', myans, 'expected:', ans, 'error:', myans-ans)

# if i use some value of dx I get my ans+eps*dx**2 where eps is some error we don't understand
# and some higher order terns = int(dx)
# if i use 2 dx I get ans+eps*(2dx)*2 = ans+4eps dx^2 = int(2dx)
# 4 * int(dx) - int(2x) = 4*ans + 4*eps*dx**2 - (ans+4*eps*dx^2) = 3*ans+.. (higher)
# so i can use ans=(4*int(dx)-int(2dx))/3

int1=myans
int2=(0.5*(y[0]+y[-1])+np.sum(y[2:-2:2]))*(2*dx)
# takes the end points and every second point
print(int1,int2)
newans=(4*int1-int2)/3 
print(newans, 'error: ', newans-ans)

# (0.5*y[0]+y[1]+y[2]...)*dx for int1
# (0.5*y[0]+y[2]+y[4]...)*2dx for int2
# (4*int1-int2)/3 = (2*y[0]+4y[1]+4y[2]..)*dx - (y[0]+2y[2]+2y[4]...)*dx 
# goes to (y[0]+4y[1]+2y[2]+4y[3]+2y[4]+...+y[-1])*dx/3

check=(y[0]+y[-1]+4*np.sum(y[1:-1:2])+2*np.sum(y[2:-2:2]))*dx/3 
print('alt answer: ', check) # we get exactly same answer and expression above

