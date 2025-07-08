from q3d.uflplus import *
import yaml

# General UFL objects

def from_director(director: list | ListTensor | Zero) -> ListTensor | Zero:
    """ vector3d -> qvector """
    return qvectorfy(qtensor_from_director(as_vector(director)))

def from_spherical_director(director: list | ListTensor | Zero) -> ListTensor:
    """ vector3d -> qvector """
    return qvectorfy(qtensor_from_director(spherical_to_cartesian(as_vector(director))))

def spherical_to_cartesian(vector_spherical: ListTensor | Zero) -> ListTensor | Zero:
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

# Specific UFL objects

def GreatCircle(*, equator_angle=pi/2, stripe_radius=pi/8):
    """ Corresponds to the Great Circle, i.e. a stripe running along the circumference of a sphere.
     
    Arguments:
        equator_angle: angle of the equator, default is pi/2 (i.e. the equator of the sphere)
        stripe_radius: radius of the stripe in terms of angle, default is pi/8

    Returns:
        A ListTensor representing the Great Circle, which is a director field that is defined as follows:
        - Outside the stripe, the director is the constant 1 and points in the r-axis (i.e. the outward pointing radial vector).
        - Inside the stripe, the director is a function of theta, which is the angle from the positive z-axis, and is given by:
          - cos(scale_factor*theta - equator_angle) along the r-axis
          - sin(scale_factor*theta - equator_angle) along the theta-axis
          - 0 along the phi-axis
          where scale_factor is pi/2 divided by the stripe_radius. This ensures that one fourth of the period of the sine and cosine functions is equal to the stripe_radius.
     """
    
    scale_factor = pi/2 / stripe_radius  # scale factor to convert stripe_radius to the range of theta

    return conditional(ge(abs(theta - equator_angle), stripe_radius), 1, 0) * from_spherical_director([1,0,0]) + conditional(lt(abs(theta - equator_angle), stripe_radius), 1, 0) * from_spherical_director([cos(scale_factor*(theta - equator_angle) + pi/2), sin(scale_factor*(theta - equator_angle) + pi/2),0])

# For the following functions, the requirement that the input is a ListTensor or Zero is no longer enforced, since it causes issues
# when defining UFL objects that are, e.g. _sums_ of ListTensors and Zeros

def add_ufl_constructors():
    def qvector_constructor(loader, node):
        def qvector(vector5d):
            # if not isinstance(vector5d, ListTensor | Zero):
            #     raise TypeError(f'Qvectors must be of type ufl.tensors.ListTensor or ufl.constantvalue.Zero, not {type(vector5d)}')
            if vector5d.ufl_shape != (5,):
                raise ValueError(f'Qvectors must have shape (5,), not {vector5d.ufl_shape}')
            return vector5d
        
        ufl_vector_object = eval(loader.construct_scalar(node))
        return qvector(ufl_vector_object)

    def vector_constructor(loader, node):
        def vector(vector3d):
            # if not isinstance(vector3d, ListTensor | Zero):
            #     raise TypeError(f'Directors must be of type ufl.tensors.ListTensor or ufl.constantvalue.Zero, not {type(vector3d)}')
            if vector3d.ufl_shape != (3,):
                raise ValueError(f'Directors must have shape (3,), not {vector3d.ufl_shape}')
            return vector3d
        
        ufl_vector_object = eval(loader.construct_scalar(node))
        return vector(ufl_vector_object)

    yaml.add_constructor('!qvector', qvector_constructor)
    yaml.add_constructor('!director', vector_constructor) # TODO: rename vector_constructor to director_constructor