from q3d.uflplus import *
import yaml

domain = AbstractDomain(3,3) # acts as placeholder domain
x0, x1, x2 = SpatialCoordinate(domain)

def from_director(director: list | ListTensor | Zero) -> ListTensor | Zero:
    """ vector3d -> qvector """
    return qvectorfy(qtensor_from_director(as_vector(director)))

def from_spherical_director(director: list | ListTensor | Zero) -> ListTensor:
    """ vector3d -> qvector """
    return qvectorfy(qtensor_from_director(spherical_to_cartesian(as_vector(director))))

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

def add_ufl_constructors():
    def qvector_constructor(loader, node):
        def qvector(vector5d: ListTensor | Zero) -> ListTensor | Zero:
            if not isinstance(vector5d, ListTensor | Zero):
                raise TypeError(f'Qvectors must be of type ufl.tensors.ListTensor or ufl.constantvalue.Zero, not {type(vector5d)}')
            if vector5d.ufl_shape != (5,):
                raise ValueError(f'Qvectors must have shape (5,), not {vector5d.ufl_shape}')
            return vector5d
        
        ufl_vector_object = eval(loader.construct_scalar(node))
        return qvector(ufl_vector_object)

    def vector_constructor(loader, node):
        def vector(vector3d: ListTensor | Zero) -> ListTensor | Zero:
            if not isinstance(vector3d, ListTensor | Zero):
                raise TypeError(f'Vectors must be of type ufl.tensors.ListTensor or ufl.constantvalue.Zero, not {type(vector3d)}')
            if vector3d.ufl_shape != (3,):
                raise ValueError(f'Vectors must have shape (3,), not {vector3d.ufl_shape}')
            return vector3d
        
        ufl_vector_object = eval(loader.construct_scalar(node))
        return vector(ufl_vector_object)

    yaml.add_constructor('!qvector', qvector_constructor)
    yaml.add_constructor('!vector', vector_constructor)