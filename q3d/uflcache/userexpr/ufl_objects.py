from q3d.ufloperatorsplus import *
import yaml

constructor_names = ('qvector_from_director', 'qvector_from_spherical_director')

def _vectorfy(tensor: ListTensor) -> ListTensor:
    from math import sqrt

    a: float = (sqrt(3.0)-3.0)/6.0
    b: float = (sqrt(3.0)+3.0)/6.0
    c: float = -sqrt(3.0)/3.0
    d: float = sqrt(2.0)/2.0

    E: list[ListTensor] = [diag(as_vector([a,b,c])), diag(as_vector([b,a,c])), as_matrix([[0,d,0],[d,0,0],[0,0,0]]), as_matrix([[0,0,d],[0,0,0],[d,0,0]]), as_matrix([[0,0,0],[0,0,d],[0,d,0]])]

    m: ListTensor = as_vector([inner(E_i,tensor) for E_i in E])
    return m

def _qtensor_from_director(director: ListTensor) -> ListTensor:
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

def _qtensor_from_spherical_director(director: ListTensor) -> ListTensor:
    from ufl import Constant, replace
    try:
        from q3d.config import constants
        S0: float = constants.S0
    except ImportError:
        S0 = 1
    
    r = Constant(domain)
    theta = Constant(domain)
    phi = Constant(domain)

    r_theta_phi = director # get r, theta, and phi coord of initial condition
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
    n = r_theta_phi[0]*e_r + r_theta_phi[1]*e_theta + r_theta_phi[2]*e_phi

    # if the director is 0, then let the Q-tensor be 0
    if inner(n,n) != 0:
        M: ListTensor = S0*(outer(n,n)/inner(n,n) - 1/3*Identity(3))
        M = conditional(ne(inner(n,n),0), M, Zero((3,3)))
    else:
        M: Zero = Zero((3,3))
    
    return M

def qvector(vector5d: list | ListTensor | Zero) -> ListTensor | Zero:
    return as_vector(vector5d)

def qvector_from_director(director: list | ListTensor | Zero) -> ListTensor | Zero:
    return _vectorfy(_qtensor_from_director(as_vector(director)))

def qvector_from_spherical_director(director: list | ListTensor | Zero) -> ListTensor:
    return _vectorfy(_qtensor_from_spherical_director(as_vector(director)))

def qvector_constructor(loader, node):
    ufl_vector_object = eval(loader.construct_scalar(node))
    return qvector(ufl_vector_object)

def add_ufl_constructors():
    yaml.add_constructor('!qvector', qvector_constructor)