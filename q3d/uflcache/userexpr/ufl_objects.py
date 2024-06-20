from q3d.ufloperatorsplus import *
import yaml

def qtensor_from_director(director: as_vector):
    try:
        from q3d.config import constants
        S0 = constants.S0
    except ImportError:
        S0 = 1
    
    n = director # get director
    n = n / sqrt(inner(n,n)) # normalize
    M = S0*(outer(n,n) - 1/3*Identity(3))
    return M

def qvector_from_director(director: as_vector):
    a = (sqrt(3.0)-3.0)/6.0
    b = (sqrt(3.0)+3.0)/6.0
    c = -sqrt(3.0)/3.0
    d = sqrt(2.0)/2.0

    E = [diag(as_vector([a,b,c])), diag(as_vector([b,a,c])), as_matrix([[0,d,0],[d,0,0],[0,0,0]]), as_matrix([[0,0,d],[0,0,0],[d,0,0]]), as_matrix([[0,0,0],[0,0,d],[0,d,0]])]

    M = qtensor_from_director(director)
    m = as_vector([inner(E_i,M) for E_i in E])
    return m

def constructor(function):
    def generic_constructor(loader, node):
        constructed_sequence = loader.construct_sequence(node)
        vector = as_vector([eval(component) for component in constructed_sequence])
        return function(vector)
    return generic_constructor

def add_ufl_constructors():
    yaml.add_constructor('!qvector_from_director', constructor(qvector_from_director))