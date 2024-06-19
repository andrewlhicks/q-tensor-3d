""" Intended to be where constructors for user-defined expressions live. """

from sympyplus import *
import yaml
from q3d.loaddump import *
from q3d.ufloperatorsplus import * # MESSES UP SYMPY SQRT FUNCTION!!!

constructors = ('SpatialVector','FromTensor','FromVector','FromDirector','FromSphericalDirector',
    'FromTensorStrongF','FromVectorStrongF','FromDirectorStrongF',
    'FromTensorStrongG','FromVectorStrongG','FromDirectorStrongG')

# Base class

class FromSympy(Matrix):
    def __new__(cls,matrix:list):
        return super().__new__(cls,matrix)
    @property
    def per_se_expr(self):
        return None
    @property
    def result(self):
        return self.per_se_expr
    def uflfy(self):
        return uflfy(self.result)
def fs_constructor(cls):
    def constructor(loader,node):
        value = loader.construct_sequence(node)
        return cls(value)
    return constructor
def fs_representer(tag):
    def representer(dumper,data):
        raise NotImplementedError('YAML dump not implemented for user expressions.')
        # data = map(str,data)
        # data = ','.join(data)
        # data = '[' + data + ']'
        # return dumper.represent_scalar(tag,data)
    return representer

# Additional base classes

class FromSympyStrongF(FromSympy):
    @property
    def result(self):
        from q3d.compute_terms import strong_F
        return vectorfy(strong_F(q_tensorfy(self.per_se_expr)))

class FromSympyStrongG(FromSympy):
    @property
    def result(self):
        from q3d.compute_terms import strong_G
        return vectorfy(strong_G(q_tensorfy(self.per_se_expr)))

# Base class for 3d vector

class SpatialVector(FromSympy):
    def __new__(cls,matrix:list):
        if len(matrix) != 3:
            raise ValueError('Vector must be length 3.')
        return super().__new__(cls,matrix)
    @property
    def per_se_expr(self):
        return self
    def __repr__(self):
        return f'FromVector {self.result}'

# Base classes for 5d vectors
class FromTensor(FromSympy):
    @property
    def per_se_expr(self):
        M = self.reshape(3,3)
        M = SymmetricTracelessMatrix(M)
        m = vectorfy(M)
        return m
    def __repr__(self):
        return f'FromTensor {self.result}'

class FromVector(FromSympy):
    def __new__(cls,matrix:list):
        if len(matrix) != 5:
            raise ValueError('FromVector must be length 5.')
        return super().__new__(cls,matrix)
    @property
    def per_se_expr(self):
        return self
    def __repr__(self):
        return f'FromVector {self.result}'

class FromDirector(FromSympy):
    @property
    def per_se_expr(self):
        try:
            from q3d.config import constants
            S0 = constants.S0
        except ImportError:
            S0 = 1
        
        n = Matrix(self)
        n = n/sqrt(n[0]**2+n[1]**2+n[2]**2)
        M = S0*(outerp(n,n) - 1/3*eye(3))
        m = vectorfy(M)
        return m
    def __repr__(self):
        return f'FromDirector {self.result}'

class FromSphericalDirector(FromSympy):
    @property
    def per_se_expr(self):
        try:
            from q3d.config import constants
            S0 = constants.S0
        except ImportError:
            S0 = 1
        
        r, theta, phi = symbols('r,theta,phi')

        r_theta_phi = Matrix(self) # get r, theta, and phi coord of initial condition
        r_theta_phi = r_theta_phi.subs(r, sqrt(x[0]**2+x[1]**2+x[2]**2))
        r_theta_phi = r_theta_phi.subs(theta, acos(x[2]/sqrt(x[0]**2+x[1]**2+x[2]**2)))
        r_theta_phi = r_theta_phi.subs(phi, sign(x[1])*acos(x[0]/sqrt(x[0]**2+x[1]**2)))

        # spherical basis vectors
        e_r = Matrix([x[0], x[1], x[2]])/sqrt(x[0]**2+x[1]**2+x[2]**2)
        e_theta = Matrix([x[0]*x[2], x[1]*x[2], -(x[0]**2+x[1]**2)])/sqrt((x[0]**2+x[1]**2+x[2]**2)*(x[0]**2+x[1]**2))
        e_phi = Matrix([-x[1], x[0], 0])/sqrt(x[0]**2+x[1]**2)

        e_theta.reshape(3,1)
        e_phi.reshape(3,1)

        # define n as a linear combination of the three basis vectors and scalars
        n = r_theta_phi[0]*e_r + r_theta_phi[1]*e_theta + r_theta_phi[2]*e_phi

        n = n/sqrt(n[0]**2+n[1]**2+n[2]**2)
        M = S0*(outerp(n,n) - 1/3*eye(3))
        m = vectorfy(M)
        return m
    def __repr__(self):
        return f'FromDirector {self.result}'

# Base classes for 5d vectors with domain strong form

class FromTensorStrongF(FromTensor,FromSympyStrongF):
    def __repr__(self):
        return f'FromTensorStrongF {self.result}'

class FromVectorStrongF(FromVector,FromSympyStrongF):
    def __repr__(self):
        return f'FromVectorStrongF {self.result}'

class FromDirectorStrongF(FromDirector,FromSympyStrongF):
    def __repr__(self):
        return f'FromDirectorStrongF {self.result}'

# Base classes for 5d vectors with boundary strong form

class FromTensorStrongG(FromTensor,FromSympyStrongG):
    def __repr__(self):
        return f'FromTensorStrongG {self.result}'

class FromVectorStrongG(FromVector,FromSympyStrongG):
    def __repr__(self):
        return f'FromVectorStrongG {self.result}'

class FromDirectorStrongG(FromDirector,FromSympyStrongG):
    def __repr__(self):
        return f'FromDirectorStrongG {self.result}'

for constructor in constructors:
    constructor_name = '!' + constructor
    constructor = eval(constructor)
    yaml.add_constructor(constructor_name,fs_constructor(constructor))
    # yaml.add_representer(constructor,fs_representer(constructor_name)) # Won't add since not implemented

### new direct UFL constructors

from ufl.tensors import ComponentTensor
from ufl.geometry import SpatialCoordinate

class FromUFL(ComponentTensor):
    pass

def qmatrix_constructor(cls):
    def constructor(loader, node):
        value = loader.construct_sequence(node)
        director = as_vector([eval(component) for component in value])
        return qvector_from_director(director)
    return constructor

yaml.add_constructor('!QmatrixFromDirector', qmatrix_constructor(FromUFL))