#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 17:26:01 2025

@author: Grace
"""

import numpy as np
import matplotlib.pyplot as plt

def lorentz(x):
    return 1/(1+x**2)

def get_weights(n):
    x = np.linspace(-1, 1, n+1)
    P = np.zeros([n+1, n+1])
    P[:,0] = 1
    P[:,1] = x
    for nn in range(1,n): # evaluate recurrance relation to get P_n
        P[:,nn+1]=((2*nn+1)*x*P[:,nn]-nn*P[:,nn-1])/(nn+1)
    Pinv = np.linalg.inv(P)
    wt = Pinv[0,:]
    # correction for total sum
    return wt/wt.sum()*n  # n+1 points for n intervals

def integrate_leg(f,a,b,dx_targ=0.1,ord=4):
    # can decide dx and order given integration limits
    # need right integer number of spacings (intervals) between starting points
    
    # do simpson's rule first (need odd number of points)
    # need 3 points to fit a quad, adding regions so need 2 more points
    npt_rough = (b-a)/dx_targ
    n_group = int(np.round(npt_rough/ord)) # get approx number of order-sized blocks
    npt = n_group*ord+1 # make sure we now have exactly correct number of points
    #print(npt) # tells us how closest number of odd number points we can get determined by spacing
    wts = get_weights(ord) # weights for legendre inetgration
    x = np.linspace(a,b,npt) # get x values
    dx = x[1] - x[0] # slightly sloppy version fo dx because might be different than average
    tot = 0
    y = f(x) # function values
    # sum weights times function values over each ord*dx interval
    for i in range(n_group):
        i_start = i*ord # each order chunk starts at some integer times the order
        tot = tot+np.sum(wts*y[i_start:i_start+ord+1]) # simpson's uses 3 poitns at a time
    return tot*dx

print('second order weights are: ', get_weights(2))
#ord = 8
#print('for order ', 2, 'weights are ', get_weights(ord)*ord)
#print('weight sum: ', np.sum(get_weights(ord)))
a = 0
b = 1
dx = 0.086
integral = integrate_leg(np.exp,a,b,dx,2)
print('')
print('integral is: ',integral)
truth = np.exp(b) - np.exp(a)
print('true value: ', truth)
print('error: ', integral - truth)

print('')

a = -1
b = 1
dx = 0.2
integral = integrate_leg(lorentz,a,b,dx,4)
print('integral is: ',integral)
truth = np.arctan(b) - np.arctan(a)
print('true value: ', truth)
print('error: ', integral - truth)





