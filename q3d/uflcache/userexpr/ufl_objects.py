from q3d.ufloperatorsplus import *
import yaml

__all__ = ('add_ufl_constructors')

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