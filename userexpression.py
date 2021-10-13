""" Intended to be where constructors for user-defined expressions live. """

from sympyplus import *
import yaml
from loaddump import *

path = 'saves/experiment/eqndata.json'
constants = load_yml(path)['constants']
S0 = constants['S0']

class UXFromTensor(Matrix):
    def __new__(cls,tensor):
        M = Matrix(tensor).reshape(3,3)
        M = SymmetricTracelessMatrix(M)
        m = vectorfy(M)
        return super().__new__(cls,m)
    def __init__(self,tensor):
        self.tensor = Matrix(tensor).reshape(3,3)
    def __repr__(self):
        return f'UXFromTensor {self.tensor}'
    def uflfy(self):
        return uflfy(self)

def uxft_constructor(loader,node):
    value = loader.construct_sequence(node)
    return UXFromTensor(value)

yaml.add_constructor('!UXFromTensor',uxft_constructor)

class UXFromVector(Matrix):
    def __new__(cls,vector):
        return super().__new__(cls,vector)
    def __repr__(self):
        return f'UXFromVector {self}'
    def uflfy(self):
        return uflfy(self)

def uxfv_constructor(loader,node):
    value = loader.construct_sequence(node)
    return UXFromVector(value)

yaml.add_constructor('!UXFromVector',uxfv_constructor)

class UXFromDirector(Matrix):
    def __new__(cls,director):
        import numpy as np
        n = Matrix(director)
        n = n/sqrt(n[0]**2+n[1]**2+n[2]**2)
        M = S0*(outerp(n,n) - 1/3*eye(3))
        m = vectorfy(M)
        return super().__new__(cls,m)
    def __init__(self,director):
        self.director = Matrix(director)
    def __repr__(self):
        return f'UXFromDirector {self.director}'
    def uflfy(self):
        return uflfy(self)

def uxfd_constructor(loader,node):
    value = loader.construct_sequence(node)
    return UXFromDirector(value)

yaml.add_constructor('!UXFromDirector',uxfd_constructor)
