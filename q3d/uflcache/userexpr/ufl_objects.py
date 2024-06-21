from q3d.ufloperatorsplus import *
import yaml

constructor_names = ('qvector', 'qvector_from_director', 'qvector_from_spherical_director')

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
    n = n / sqrt(inner(n,n)) # normalize
    M: ListTensor = S0*(outer(n,n) - 1/3*Identity(3)) # 
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
    r_theta_phi = replace(r_theta_phi, {r: sqrt(x0**2+x1**2+x2**2)})
    r_theta_phi = replace(r_theta_phi, {theta: acos(x2/sqrt(x0**2+x1**2+x2**2))})
    r_theta_phi = replace(r_theta_phi, {phi: sign(x1)*acos(x0/sqrt(x0**2+x1**2))})

    # spherical basis vectors
    e_r = as_vector([x0, x1, x2])/sqrt(x0**2+x1**2+x2**2)
    e_theta = as_vector([x0*x2, x1*x2, -(x0**2+x1**2)])/sqrt((x0**2+x1**2+x2**2)*(x0**2+x1**2))
    e_phi = as_vector([-x1, x0, 0])/sqrt(x0**2+x1**2)

    # define n as a linear combination of the three basis vectors and scalars
    n = r_theta_phi[0]*e_r + r_theta_phi[1]*e_theta + r_theta_phi[2]*e_phi
    n = n / sqrt(inner(n,n)) # normalize
    M: ListTensor = S0*(outer(n,n) - 1/3*Identity(3)) # 
    return M

def qvector(vector5d: ListTensor | Zero) -> ListTensor | Zero:
    return vector5d

def qvector_from_director(director: ListTensor | Zero) -> ListTensor | Zero:
    return _vectorfy(_qtensor_from_director(director))

def qvector_from_spherical_director(director: ListTensor) -> ListTensor:
    return _vectorfy(_qtensor_from_spherical_director(director))

def constructor(func):
    def generic_constructor(loader, node):
        constructed_sequence = loader.construct_sequence(node)
        vector: ListTensor = as_vector([eval(component) for component in constructed_sequence])
        return func(vector)
    return generic_constructor

def add_ufl_constructors():
    for constructor_name in constructor_names:
        yaml.add_constructor(f'!{constructor_name}', constructor(eval(constructor_name)))