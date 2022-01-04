""" Intended to be where constructors for user-defined expressions live. """

from sympyplus import *
import yaml
from loaddump import *

class FromTensor(Matrix):
    def __new__(cls,tensor):
        M = Matrix(tensor).reshape(3,3)
        M = SymmetricTracelessMatrix(M)
        m = vectorfy(M)
        return super().__new__(cls,m)
    def __init__(self,tensor):
        self.tensor = Matrix(tensor).reshape(3,3)
    def __repr__(self):
        return f'FromTensor {self.tensor}'
    def uflfy(self):
        return uflfy(self)

def ft_constructor(loader,node):
    value = loader.construct_sequence(node)
    return FromTensor(value)

yaml.add_constructor('!FromTensor',ft_constructor)

class FromVector(Matrix):
    def __new__(cls,vector):
        return super().__new__(cls,vector)
    def __repr__(self):
        return f'FromVector {self}'
    def uflfy(self):
        return uflfy(self)

def fv_constructor(loader,node):
    value = loader.construct_sequence(node)
    return FromVector(value)

yaml.add_constructor('!FromVector',fv_constructor)

class FromDirector(Matrix):
    def __new__(cls,director):
        from config import constants as c
        n = Matrix(director)
        n = n/sqrt(n[0]**2+n[1]**2+n[2]**2)
        M = c.S0*(outerp(n,n) - 1/3*eye(3))
        m = vectorfy(M)
        return super().__new__(cls,m)
    def __init__(self,director):
        self.director = Matrix(director)
    def __repr__(self):
        return f'FromDirector {self.director}'
    def uflfy(self):
        return uflfy(self)

def fd_constructor(loader,node):
    value = loader.construct_sequence(node)
    return FromDirector(value)

yaml.add_constructor('!FromDirector',fd_constructor)
