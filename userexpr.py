""" Intended to be where constructors for user-defined expressions live. """

from sympyplus import *
import yaml
from loaddump import *

class FromSympy(Matrix):
    def __new__(cls,matrix:list):
        return super().__new__(cls,matrix)
    @property
    def result(self):
        return None
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

class FromTensor(FromSympy):
    @property
    def result(self):
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
    def result(self):
        return self
    def __repr__(self):
        return f'FromVector {self.result}'
yaml.add_constructor('!FromVector',fs_constructor(FromVector))
yaml.add_representer(FromVector,fs_representer('!FromVector'))

class FromDirector(FromSympy):
    @property
    def result(self):
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