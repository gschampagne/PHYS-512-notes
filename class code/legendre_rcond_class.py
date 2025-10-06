#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  5 13:14:52 2025

@author: Grace
"""

import numpy as np

x = np.linspace(-1,1,1001)
order = [5,10,15,20,30,50,75,100]
for ord in order:
    A = np.polynomial.legendre.legvander(x,ord)
    u,s,v = np.linalg.svd(A)
    print('for order ', ord, ' condition is ', s.max()/s.min())