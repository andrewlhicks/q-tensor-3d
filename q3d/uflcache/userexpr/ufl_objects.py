from firedrake import BoxMesh
from q3d.uflplus import *
from ufl import as_vector as _as_vector
import yaml

domain = BoxMesh(1, 1, 1, 1, 1, 1) # acts as completely arbitrary placeholder domain (AbstractDomain fails with spherical coords)
x0, x1, x2 = SpatialCoordinate(domain)

r = Constant(domain)
theta = Constant(domain)
phi = Constant(domain)

def as_vector(vector: Expr) -> Expr:
    """ Replaces the UFL as_vector function only in this module """

    vector = _as_vector(vector)

    vector = replace(vector, {r: sqrt(x0**2+x1**2+x2**2)})
    vector = replace(vector, {theta: acos(x2/sqrt(x0**2+x1**2+x2**2))})
    vector = replace(vector, {phi: sign(x1)*acos(x0/sqrt(x0**2+x1**2))})

    return vector

def from_director(director: Expr) -> Expr:
    """ vector3d -> qvector """
    return qvectorfy(qtensor_from_director(as_vector(director)))

def from_spherical_director(director: Expr) -> ListTensor:
    """ vector3d -> qvector """
    return qvectorfy(qtensor_from_director(spherical_to_cartesian(as_vector(director))))

def spherical_to_cartesian(vector_spherical: Expr) -> Expr:
    """ vector3d -> vector3d """

    if vector_spherical.ufl_shape != (3,):
        raise ValueError(f'Spherical vector must have shape (3,), not {vector_spherical.ufl_shape}')

    # spherical basis vectors
    e_r = as_vector([x0, x1, x2])/sqrt(x0**2+x1**2+x2**2)
    e_theta = as_vector([x0*x2, x1*x2, -(x0**2+x1**2)])/sqrt((x0**2+x1**2+x2**2)*(x0**2+x1**2))
    e_phi = as_vector([-x1, x0, 0])/sqrt(x0**2+x1**2)

    # remove anything undefined
    e_r = conditional(ne(x0**2+x1**2+x2**2,0), e_r, as_vector([0,0,0]))
    e_theta = conditional(ne(x0**2+x1**2,0), e_theta, as_vector([0,0,0]))
    e_phi = conditional(ne(x0**2+x1**2,0), e_phi, as_vector([0,0,0]))

    # define n as a linear combination of the three basis vectors and scalars
    vector_cartesian = vector_spherical[0]*e_r + vector_spherical[1]*e_theta + vector_spherical[2]*e_phi

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
        def qvector(vector5d: Expr) -> Expr:
            if not isinstance(vector5d, Expr):
                raise TypeError(f'Qvectors must be of type ufl.tensors.ListTensor or ufl.constantvalue.Zero, not {type(vector5d)}')
            if vector5d.ufl_shape != (5,):
                raise ValueError(f'Qvectors must have shape (5,), not {vector5d.ufl_shape}')
            return vector5d
        
        ufl_vector_object = eval(loader.construct_scalar(node))
        return qvector(ufl_vector_object)

    def vector_constructor(loader, node):
        def vector(vector3d: Expr) -> Expr:
            if not isinstance(vector3d, Expr):
                raise TypeError(f'Vectors must be of type ufl.tensors.ListTensor or ufl.constantvalue.Zero, not {type(vector3d)}')
            if vector3d.ufl_shape != (3,):
                raise ValueError(f'Vectors must have shape (3,), not {vector3d.ufl_shape}')
            return vector3d
        
        ufl_vector_object = eval(loader.construct_scalar(node))
        return vector(ufl_vector_object)

    yaml.add_constructor('!qvector', qvector_constructor)
    yaml.add_constructor('!vector', vector_constructor)