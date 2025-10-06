#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  5 13:06:45 2025

@author: Grace
"""

import numpy as np
import matplotlib.pyplot as plt
plt.ion()

x = np.linspace(-1,1,1001)
y22 = x**22
y24 = x**24 
y26 = x**26
yy24 = (y22 + y26) / 2
plt.clf()
plt.plot(x,y24)
plt.plot(x,yy24)
plt.show()