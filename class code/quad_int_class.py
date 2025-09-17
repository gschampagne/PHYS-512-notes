#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 15:31:50 2025

@author: Grace
"""

import numpy as np

# we have f(x=-1)=yl, f(x=0)=y0, f(x=1)=yr
# we want quadratic ax^2 + bx + c that goes through these points
# c=y0
# do we care about b? no since between 1 and -1
# yl = a - b + c
# yr = a + b + c
# y0 = c

# yl + yr = 2a + 2c, yl -2y0 + yr = 2a, a = (yl + yr - 2y0)/2
# area is the integral from -1 to 1
# integral is cx + ax^3/3 from -1 to 1
# gives 2c + 2/3a = 2y0 + 2/3 (yl + yr - 2y0)/2 = 1/3 (6y0 + yl + yr - 2y0)
# = 1/3(yl + 4y0 + yr)
# a bunch of y values: first interval is 1/3 (y[0] + 4y[1] + y[2]) +
# + 1/3(y[2] + 4y[3] + y[4]) + ....
# 1/3 (y[0] + 4y[1] + 2y[2] + 4y[3] + 2y[4] + ... + y[-1])
# same alternating pattern of 1 on the ends and 2, 4, s for middle coeff
# exactly the same as the error cancellation we calculated for the linear inetgration
# called simpson's rule and is a useful integrator

def f(x):
    return 1/x


def simpson_raw(x,y):
    dx = x[1] - x[0]
    return dx*(y[0] + y[-1] + 4*np.sum(y[1:-1:2])+2*np.sum(y[2:-2:2]))/3 

def simpson(f,a,p,npt):
    x=np.linspace(a,b,npt)
    y=f(x)
    return simpson_raw(x,y)

a=1
b=2

print('1/x from simpson is', simpson(f,1,2,41), 'true answer is', np.log(b)-np.log(a))
# need odd number of points for simpson's rule
# if double number of points get factor of 16 more accuracy
# did not have to interpolate directly sicne already can integrate

