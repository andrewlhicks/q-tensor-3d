from firedrake import BoxMesh
from ufl import SpatialCoordinate, Constant, replace, acos, sqrt, sign

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