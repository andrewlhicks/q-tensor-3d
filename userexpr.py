""" Intended to be where constructors for user-defined expressions live. """

from sympyplus import *
import yaml
from loaddump import *

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
    def strong_form_domain(self,vector):
        from compute_terms import strong_F
        return vectorfy(strong_F(q_tensorfy(vector)))
    def strong_form_boundary(self,vector):
        from compute_terms import strong_G
        return vectorfy(strong_G(q_tensorfy(vector)))
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
yaml.add_constructor('!SpatialVector',fs_constructor(SpatialVector))
yaml.add_representer(SpatialVector,fs_representer('!SpatialVector'))

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
yaml.add_constructor('!FromTensor',fs_constructor(FromTensor))
yaml.add_representer(FromTensor,fs_representer('!FromTensor'))

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
yaml.add_constructor('!FromVector',fs_constructor(FromVector))
yaml.add_representer(FromVector,fs_representer('!FromVector'))

class FromDirector(FromSympy):
    @property
    def per_se_expr(self):
        from config import constants as c
        n = Matrix(self)
        n = n/sqrt(n[0]**2+n[1]**2+n[2]**2)
        M = c.S0*(outerp(n,n) - 1/3*eye(3))
        m = vectorfy(M)
        return m
    def __repr__(self):
        return f'FromDirector {self.result}'
yaml.add_constructor('!FromDirector',fs_constructor(FromDirector))
yaml.add_representer(FromDirector,fs_representer('!FromDirector'))

# Base classes for 5d vectors with domain strong form

class FromTensorStrongF(FromTensor):
    @property
    def result(self):
        return self.strong_form_domain(self.per_se_expr)
    def __repr__(self):
        return f'FromTensorStrongF {self.result}'
yaml.add_constructor('!FromTensorStrongF',fs_constructor(FromTensorStrongF))
yaml.add_representer(FromTensorStrongF,fs_representer('!FromTensorStrongF'))

class FromVectorStrongF(FromVector):
    @property
    def result(self):
        return self.strong_form_domain(self.per_se_expr)
    def __repr__(self):
        return f'FromVectorStrongF {self.result}'
yaml.add_constructor('!FromVectorStrongF',fs_constructor(FromVectorStrongF))
yaml.add_representer(FromVectorStrongF,fs_representer('!FromVectorStrongF'))

class FromDirectorStrongF(FromDirector):
    @property
    def result(self):
        return self.strong_form_domain(self.per_se_expr)
    def __repr__(self):
        return f'FromDirectorStrongF {self.result}'
yaml.add_constructor('!FromDirectorStrongF',fs_constructor(FromDirectorStrongF))
yaml.add_representer(FromDirectorStrongF,fs_representer('!FromDirectorStrongF'))

# Base classes for 5d vectors with boundary strong form

class FromTensorStrongG(FromTensor):
    @property
    def result(self):
        return self.strong_form_boundary(self.per_se_expr)
    def __repr__(self):
        return f'FromTensorStrongG {self.result}'
yaml.add_constructor('!FromTensorStrongG',fs_constructor(FromTensorStrongG))
yaml.add_representer(FromTensorStrongG,fs_representer('!FromTensorStrongG'))

class FromVectorStrongG(FromVector):
    @property
    def result(self):
        return self.strong_form_boundary(self.per_se_expr)
    def __repr__(self):
        return f'FromVectorStrongG {self.result}'
yaml.add_constructor('!FromVectorStrongG',fs_constructor(FromVectorStrongG))
yaml.add_representer(FromVectorStrongG,fs_representer('!FromVectorStrongG'))

class FromDirectorStrongG(FromDirector):
    @property
    def result(self):
        return self.strong_form_boundary(self.per_se_expr)
    def __repr__(self):
        return f'FromDirectorStrongG {self.result}'
yaml.add_constructor('!FromDirectorStrongG',fs_constructor(FromDirectorStrongG))
yaml.add_representer(FromDirectorStrongG,fs_representer('!FromDirectorStrongG'))