# Polynomial Interpolation/ Integration (class slides posted)

### Interpolation


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
