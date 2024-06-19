from ufl.operators import *
from ufl import Identity, AbstractDomain

domain = AbstractDomain(3,3)
x0, x1, x2 = SpatialCoordinate(domain)

def smooth_transition(x, a, b):
    """ A smooth transition function in the variable x wih the
    value of 0 at x=a and 1 at x=b """
    # smooth function f defined by the exponential
    def f(x):
        return conditional(gt(x,0), exp(-1/x), 0)
    # smooth function g such that g(x)=0 for x<=0 and g(x)=1 for x>=1
    def g(x):
        return f(x)/(f(x)+f(1-x))
    
    return g((x-a)/(b-a))

def qtensor_from_director(director: as_vector):
    try:
        from q3d.config import constants
        S0 = constants.S0
    except ImportError:
        S0 = 1
    
    n = director # get director
    n = n / sqrt(inner(n,n)) # normalize
    M = S0*(outer(n,n) - 1/3*Identity(3))
    return M

def qvector_from_director(director: as_vector):
    a = (sqrt(3.0)-3.0)/6.0
    b = (sqrt(3.0)+3.0)/6.0
    c = -sqrt(3.0)/3.0
    d = sqrt(2.0)/2.0

    E = [diag(as_vector([a,b,c])), diag(as_vector([b,a,c])), as_matrix([[0,d,0],[d,0,0],[0,0,0]]), as_matrix([[0,0,d],[0,0,0],[d,0,0]]), as_matrix([[0,0,0],[0,0,d],[0,d,0]])]

    M = qtensor_from_director(director)
    m = as_vector([inner(E_i,M) for E_i in E])
    return m