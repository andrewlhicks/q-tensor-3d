from q3d.ufloperatorsplus import *
import yaml

def from_director(director: list | ListTensor | Zero) -> ListTensor | Zero:
    return vectorfy(qtensor_from_director(as_vector(director)))

def from_spherical_director(director: list | ListTensor | Zero) -> ListTensor:
    return vectorfy(qtensor_from_director(spherical_to_cartesian(as_vector(director))))

def add_ufl_constructors():
    def qvector_constructor(loader, node):
        def qvector(vector5d: list | ListTensor | Zero) -> ListTensor | Zero:
            q = as_vector(vector5d)
            if q.ufl_shape != (5,):
                raise ValueError(f'Qvectors must have shape (5,), not {q.ufl_shape}')
            return q
        
        ufl_vector_object = eval(loader.construct_scalar(node))
        return qvector(ufl_vector_object)

    yaml.add_constructor('!qvector', qvector_constructor)