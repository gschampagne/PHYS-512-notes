#### wrote code for the equation above 

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


# Integration


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

