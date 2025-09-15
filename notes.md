#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 08:31:43 2025

@author: Grace
"""

#### Lecture 1 - Aug ?

# Error in the Derivative


using classical definition of derivative: f'(x) = (f(x+dx)-f(x))/dx

then taylor expand both f(x+dx) and f(x)
substitute back into equation and get: f'(x) = f'(x) + dx^2 f''/2 + ...

but now try: f'(x) = (f(x+dx)-f(x-dx))/(2*dx) 
again taylor expand and sub back in and now the dx^2 terms dissapear due to sign

f'(x) = f'(x) + dx^2 f'''/6 + ...

error in taylor series decreases by factor of 100 - this answer is way more accurate
*have to worry about round off error because now is dominating
-> want to make dx larger

so for single precision would want like 0.1 and for double would be like 10e5

summary: free to pick dx whatever we want and so to get more accurate we choose 
        dx to be the cubed root fo the third derivative (equalling sides of the 
        equation) rather than the root of the second derivative

- see code showing this called deriv_class.py
- most accurate we can get is 10e-10


#### Lecture 2 - Sept 2

# Polynomial Interpolation/ Integration (class slides posted)

## Interpolation


if have a function given at certain points and asked to find the function 
value y at some other x value
- interpolation is when the new value is between already given function values
- extrapolation is when new value is outside of given function values

solutions:
- laziest way is to choose closest neighbour
- draw straight lines between points (also lazy)
- piecewise quadratic (not good usually)
**we want a smooth function a.k.a. something that can be defined by a taylor series

taylor series reminder
f(x) = 

- cubic (needs 4 points) use edge points and so thye match up on center point
    and is smooth
        -> see code from class cubic_attempt_class_py but indexing didn't work
        - cubic works well!
        - what about higher order?

I want to write a polynomial that goes through n points
- can write a polynomial that is 0 at specified number of points
    - ex. (x0 - x1)(x0 - x2)(..)(x0 - xn-1) which is 1 at x0 and 0 at all other values
say we have function y(x) = y0@x0, y1@x1, ... = y0P0 + y1P1 + y2P2 + ...
- could be defined by sum of polynomials P which go through zero except for point
    that they 'own' and so they wiggle a lot
    - see images in lecture slides

note that which many high order polynomials there can be large effects of small
errors and so not well behaved

checkin: cubic worked ok and was continuous. did nto need fit for each 
interpolated point, could store every interval and coeff
- no need evenly spaced given function values


### Cubic Splines


our cubic fit, was it continuous? we did not check :( 

splines: forces function and first n order derivatievs to be continuous 
        (usually 1st and 2nd)
        - so second deriv needs to match right and left neighbours (usually set to 0)
        - cubic most common to be a spline
        - see code cubic_spline_class.py
        - really fancy is Bplines found in scipy.interpolate


### Rational Functions

poles review (omfg complex analysis)

ex. function: f(x) 1/ (2 + x^2) which is not analytic at x = i
- taylor series work in the complex plane (since function is analyic)
- interpolation with spline will fail at some point

f(x) = P(x)/(1 + qq(x))

if i use a polynomial and extrapolate to long distances then they will move and 
change very steeply when out of the range we fitted them in

basically have will have 

(P0 + P1 x + P2 x^2 + .. ) / ( 1 + qq1 x + qq2 x^2 + ..) = y(x)
rearrange to 
P0 + P1 x + P2 x^2 ... - ( 1 + qq1 x + qq2 x^2 + ..) = y(x)

now have a matrix and can solve for all unknowns 
(if have that same amount of function values)




#### Lecture 3 - Sept 4

wrote code for the equation above 

f(x) = P(x)/(1 + qq(x))

for a gaussian
- code called ratfit_class.py

in code we try using n=2 and m=4 but get a splot where it fits fine but there
is a spike exactly at x=0
- if change to n=3 it works well (basically perfect polynomial fit to the points)
- only rational functions can fit these points to the gaussian so well 
    (even outside the range of points)
- n=4 and n=6 also has spike

gaussian is an even function
- n=3 is the qudratic so it is even so the fit works
- n=4 is an odd polynomial so the extra odd term does not do good

if multiply top and bottom by same thing (polynomial factor)  that causes the 
spike since the inverse blows up due to an undetermined term
- but, if have something that doesn't work like n=2 and m=6 can change this one line:
    fitp=np.linalg.inv(A)@y to fitp=np.linalg.pinv(A)@y to fix

what is pinv?
- normal when inverse works it makes eigenvalues go to 1/eigenvalues
- pseudo inverse will set 1/0=0 instead of infinity so does not break

summary: rational functions often behave better outside the region of specified values
        though they do sometimes have singularies and some issues that polys don't

note: rational functions are much less supported than polyfits


## Integration


Interpolation and integration are closely couples

If i'm trying to numerically inetrgate between a set of values I want to be able to 
interpolate what happens in that interval then integrate that
- fairly easy with polynomials
- usually end up with a set of coeff times function values wehre coeff set by scheme
    you've chosen to do the integration
- can think of this as finding the 'average' value in a region, based on some
    interpolation scheme


### Integration with Linear

interpolating by drawing a straight line between points of a fucntion
- the average value would be the average of the endpoints of each region
- if spaced by dx, the area of the interval would be average*dx
    - so area/dx = 0.5(y0+y1) + 0.5(y1+y2) + ... + 0.5(yn-2+yn-1)
                = 0.5(y0 + yn-1) + sum(yi), i=1,2...n-2
                = the ends plus the middles

wrote a code to determine if this is true called interp_linear_class.py
- got pretty well (error less than a thousandth)
- how would we expect error to scale by number of points?
    - think about taylor expansion, the linear expansion gets first order term right,
        second order term survives, get third order that when integrate turns into second 
        order answer and so would expect error in tengral to go like 1/number of points^2
    - error scales like 1/(number of points^2)
    - if increase number of points by factor of two then error should go down by a factor of 4

can we use this scheme to figure out a better answer?
- if i use some value of dx I get my ans+eps*dx**2 where eps is some error we don't understand
    and some higher order terns = int(dx)
- if i use 2 dx I get ans+eps*(2dx)*2 = ans+4eps dx^2 = int(2dx)
- 4 * int(dx) - int(2x) = 4*ans + 4*eps*dx**2 - (ans+4*eps*dx^2) = 3*ans+.. (higher)
    - so i can use ans=(4*int(dx)-int(2dx))/3

wrote a code to simulate this called interp_linear_combo.py
by using this we get a new smaller error of 3.05e-6
- got a much more accurate answer with not much extra work
- error would divide by 16 in doubled amount of points
    - leading term in taylor series is now the 4th deg term
- by being more clever in how to use functions can get much smaller error


### Integration with quadratic with 3 points

see code called quad_int_class.py
- introduce simpson's rule


# Legendre Polynomials 

we want to do even better but will need legendre polynomials since regular polynomials will 
not work

we generate legendre polynomials through a recurrance relationship
- where the previous terms are used to determine the following term in the polynomial
- (n+1)P_(n+1) = (2n+1)xP_n - nP_(n-1), P_0=1, P_1=x
- this is a set of polynomials each term increasing in order

Lengendre polynomilas are:
- bound on -1 to 1 
- integral of any two polynomials are zero if n does not equal m
    - int(P_n * P_m)= delta_n,m

visualized legendre polynomials on legendre.py

want to know inetgral of P_n for any n between 1 to 1
- for P_0 the integral is 2 since P_0 is 1
- for any other legendre polynomial the integral will be zero (can see this in the code visualization)

therefore, I can fit legendre poly to a set of data and the netgral will be the P_0 coeff!

legendre fitting: y_i = sum(c_j P_j(x_i))
- if figure out c_j then we good
- this is a matrix equation -> y=Pc and P is square if we have as many polynomials as we have points
    meaning we could use c=P^-1 y to solve for c, get c0 and we have integral
- but c_0 is just sum(P^(-1)_0,k y_k) so can just take first column of P^(-1) to get weights
- now can integrate to any order! (making sure i have suitable yi for the chosen order)

see code legendre_weights.py
- can see it essentially matches simpson's rule (up to order 2)



