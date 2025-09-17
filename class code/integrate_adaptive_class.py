#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 19:22:30 2025

@author: Grace
"""

import numpy as np

def lorentz(x):
    return 1/(1+x**2)

def gauss_off(x):
    return 1+np.exp(-0.5*(x/2)**2)

def simp_step(f,a,b):
    x = np.linspace(a,b,5)
    y = f(x)
    dx = x[1] - x[0]
    ans1 = (y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])*dx/3
    ans2 = (y[0]+4*y[2]+y[4])*(2*dx)/3
    return ans1, np.abs(ans1-ans2) # return function, error estimate

# now use recursion
def simp_adaptive(f,a,b,tol):
    #print('integrating from', a, 'to', b)
    val,err = simp_step(f,a,b)
    if err<tol:
        return val
    else:
        mid = (a+b)/2
        return simp_adaptive(f,a,mid,tol/2)+simp_adaptive(f,mid,b,tol/2)

a = -20
b = 15
tol = 1e-6
ans = simp_adaptive(lorentz, a, b, tol)
print('adaptive intgeral: ', ans, 'error: ',ans-(np.arctan(b)-np.arctan(a)))

ans = simp_adaptive(gauss_off, a, b, tol)
print('adaptive intgeral: ', ans, 'error: ',ans-(np.arctan(b)-np.arctan(a)))

