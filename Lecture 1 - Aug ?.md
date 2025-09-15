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