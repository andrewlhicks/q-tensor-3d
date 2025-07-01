from ufl import conditional, exp, gt

def smooth_transition(x, *, I: list[int | float]):
    """ A smooth transition function in the variable x on the
    interval I = [a,b], with the value of 0 at x=a and 1 at x=b """

    a = I[0]
    b = I[1]

    # smooth function f defined by the exponential
    def f(x):
        return conditional(gt(x,0), exp(-1/x), 0)
    # smooth function g such that g(x)=0 for x<=0 and g(x)=1 for x>=1
    def g(x):
        return f(x)/(f(x)+f(1-x))
    
    return g((x-a)/(b-a))

# def reverse_smooth_transition(x, *, I: list[int | float]):
#     """ A smooth transition function in the variable x on the
#     interval I = [a,b], with the value of 1 at x=a and 0 at x=b """

#     a = I[0]
#     b = I[1]

#     # smooth function f defined by the exponential
#     def f(x):
#         return conditional(gt(x,0), exp(-1/x), 0)
#     # smooth function g such that g(x)=0 for x<=0 and g(x)=1 for x>=1
#     def g(x):
#         return f(1-x)/(f(1-x)+f(x))
    
#     return g((x-a)/(b-a))

def reverse_smooth_transition(x, *, I: list[int | float]):
    return 1 - smooth_transition(x, I=I)

def gaussian(x, *, mu=0, sigma=1):
    """ Gaussian function with mean mu and standard deviation sigma """
    return exp(-(x - mu)**2 / (2 * sigma**2))