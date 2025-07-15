from firedrake import BoxMesh
from ufl import SpatialCoordinate, Constant, replace, acos, sqrt, sign
import functools

domain = BoxMesh(1, 1, 1, 1, 1, 1) # acts as completely arbitrary placeholder domain (AbstractDomain fails with spherical coords)
x0, x1, x2 = SpatialCoordinate(domain)

r = Constant(domain)
theta = Constant(domain)
phi = Constant(domain)

def process_spherical_coordinates(expr):
    """ Process spherical coordinates by converting to raw Cartesian form. """
    expr = replace(expr, {r: sqrt(x0**2+x1**2+x2**2),
                          theta: acos(x2/sqrt(x0**2+x1**2+x2**2)),
                          phi: sign(x1)*acos(x0/sqrt(x0**2+x1**2))})
    return expr

def spherical_compatible(func):
    """ Decorator that makes a function spherical-compatible by replacing spherical coordinates with Cartesian ones. """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        args_0 = process_spherical_coordinates(args[0])
        return func(args_0, *args[1:], **kwargs)
    return wrapper