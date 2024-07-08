from ufl import as_vector, as_matrix, diag, inner
from ufl.core.expr import Expr
from math import sqrt

a: float = (sqrt(3.0)-3.0)/6.0
b: float = (sqrt(3.0)+3.0)/6.0
c: float = -sqrt(3.0)/3.0
d: float = sqrt(2.0)/2.0

E: list = [diag(as_vector([a,b,c])), diag(as_vector([b,a,c])), as_matrix([[0,d,0],[d,0,0],[0,0,0]]), as_matrix([[0,0,d],[0,0,0],[d,0,0]]), as_matrix([[0,0,0],[0,0,d],[0,d,0]])]

def qvectorfy(tensor: Expr) -> Expr:
    """ qtensor -> qvector """
    return as_vector([inner(E_i,tensor) for E_i in E])

def qtensorfy(vector: Expr) -> Expr:
    """ qvector -> qtensor """
    return sum([vector_i*E_i for vector_i, E_i in zip(vector,E)])
