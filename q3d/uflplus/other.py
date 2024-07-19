from ufl import conditional, exp, gt

def smooth_transition(x, *, I: list):
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