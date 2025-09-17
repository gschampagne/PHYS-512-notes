#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 18:54:22 2025

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
    npt_rough = (b-a)/dx_targ
    n_group = int(np.round(npt_rough/ord))
    return integrate_leg_worker(f,a,b,n_group,ord)

def integrate_leg_worker(f,a,b,n_group,ord):
    wts = get_weights(ord)
    npt = n_group*ord+1
    x = np.linspace(a,b,npt)
    y=f(x)
    dx = x[1] - x[0]
    tot = 0
    for i in range(n_group):
        i_start = i*ord
        tot = tot+np.sum(wts*y[i_start:i_start+ord+1])
    return tot*dx

def integrate_leg_werr(f,a,b,dx_targ=0.1,ord=4):
    # first, make sure we have consistent number of points given our order
    npt_rough = (b-a)/dx_targ
    n_group = int(np.round(npt_rough/ord)) # get approx number of order sized blocks
    # make sure have even number of groups
    if n_group%2==1:
        n_group=n_group+1
    ans_fine = integrate_leg_worker(f,a,b,n_group,ord)
    ans_coarse = integrate_leg_worker(f, a, b, n_group//2, ord) # half as many points
    return ans_fine, np.abs(ans_fine - ans_coarse)
    






