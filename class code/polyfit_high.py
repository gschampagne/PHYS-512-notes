import numpy as np
import matplotlib.pyplot as plt

t = np.arange(-5,5,0.1)
x_true = t**3 - 0.5*t**2
x = x_true + 10*np.random.randn(t.size)

npoly = 23 # fit 4th order polynomial
ndata = t.size
A = np.zeros([ndata,npoly])
A[:,0] = 1.0
for i in range(1,npoly):
    A[:,i] = A[:,i-1]*t

# let's ignore noise for now
# new eq are m = (A^TA)^(-1)(A^Td)
A = np.matrix(A)
d = np.matrix(x).transpose()
lhs = A.transpose()*A
rhs = A.transpose()*d
fitp = np.linalg.inv(lhs)*rhs
pred = A*fitp

plt.clf();plt.plot(t,x,'*');plt.plot(t,pred,'r')
plt.draw()