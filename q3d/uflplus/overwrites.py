from ufl.operators import ListTensor, Zero
from ufl import as_vector as _as_vector
from ufl import conditional as _conditional
from q3d.uflplus.domain import process_spherical_coordinates, spherical_compatible

# @spherical_compatible
# CANNOT use the @spherical_compatible decorator here because _as_vector must be called before spherical_compatible is applied
def as_vector(vector: list | ListTensor | Zero) -> ListTensor | Zero:
    """ Replaces the UFL as_vector function only in this module """

    vector = _as_vector(vector)
    vector = process_spherical_coordinates(vector)

    return vector

@spherical_compatible
def conditional(condition, true_value, false_value):
    """ Replaces the UFL conditional function only in this module """
    return _conditional(condition, true_value, false_value)