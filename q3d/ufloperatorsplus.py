from ufl.operators import *
from ufl.classes import *
from ufl import Identity, AbstractDomain
from firedrake.ufl_expr import *

domain = AbstractDomain(3,3)
x0, x1, x2 = SpatialCoordinate(domain)

def from_director(director: list | ListTensor | Zero) -> ListTensor | Zero:
    """ vector3d -> qvector """
    return vectorfy(qtensor_from_director(as_vector(director)))

def from_spherical_director(director: list | ListTensor | Zero) -> ListTensor:
    """ vector3d -> qvector """
    return vectorfy(qtensor_from_director(spherical_to_cartesian(as_vector(director))))

def spherical_to_cartesian(vector_spherical: ListTensor | Zero) -> ListTensor | Zero:
    """ vector3d -> vector3d """
    from ufl import Constant, replace

    if vector_spherical.ufl_shape != (3,):
        raise ValueError(f'Spherical vector must have shape (3,), not {vector_spherical.ufl_shape}')

    r = Constant(domain)
    theta = Constant(domain)
    phi = Constant(domain)

    r_theta_phi = vector_spherical # get r, theta, and phi coord of initial condition
    # r_theta_phi = replace(r_theta_phi, {r: sqrt(x0**2+x1**2+x2**2)})
    # r_theta_phi = replace(r_theta_phi, {theta: acos(x2/sqrt(x0**2+x1**2+x2**2))})
    # r_theta_phi = replace(r_theta_phi, {phi: sign(x1)*acos(x0/sqrt(x0**2+x1**2))})

    # spherical basis vectors
    e_r = as_vector([x0, x1, x2])/sqrt(x0**2+x1**2+x2**2)
    e_theta = as_vector([x0*x2, x1*x2, -(x0**2+x1**2)])/sqrt((x0**2+x1**2+x2**2)*(x0**2+x1**2))
    e_phi = as_vector([-x1, x0, 0])/sqrt(x0**2+x1**2)

    # remove anything undefined
    e_r = conditional(ne(x0**2+x1**2+x2**2,0), e_r, as_vector([0,0,0]))
    e_theta = conditional(ne(x0**2+x1**2,0), e_theta, as_vector([0,0,0]))
    e_phi = conditional(ne(x0**2+x1**2,0), e_phi, as_vector([0,0,0]))

    # define n as a linear combination of the three basis vectors and scalars
    vector_cartesian = r_theta_phi[0]*e_r + r_theta_phi[1]*e_theta + r_theta_phi[2]*e_phi

    return vector_cartesian

def qtensor_from_director(director: ListTensor) -> ListTensor:
    """ vector3d -> qtensor """
    try:
        from q3d.config import constants
        S0: float = constants.S0
    except ImportError:
        S0 = 1
    
    n: ListTensor = director # get director

    # if the director is 0, then let the Q-tensor be 0
    if inner(n,n) != 0:
        M: ListTensor = S0*(outer(n,n)/inner(n,n) - 1/3*Identity(3))
        M = conditional(ne(inner(n,n),0), M, Zero((3,3)))
    else:
        M: Zero = Zero((3,3))
    
    return M

def vectorfy(tensor: ListTensor) -> ListTensor:
    """ qtensor -> qvector """
    from math import sqrt

    a: float = (sqrt(3.0)-3.0)/6.0
    b: float = (sqrt(3.0)+3.0)/6.0
    c: float = -sqrt(3.0)/3.0
    d: float = sqrt(2.0)/2.0

    E: list[ListTensor] = [diag(as_vector([a,b,c])), diag(as_vector([b,a,c])), as_matrix([[0,d,0],[d,0,0],[0,0,0]]), as_matrix([[0,0,d],[0,0,0],[d,0,0]]), as_matrix([[0,0,0],[0,0,d],[0,d,0]])]

    m: ListTensor = as_vector([inner(E_i,tensor) for E_i in E])
    return m

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