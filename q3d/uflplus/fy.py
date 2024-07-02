from ufl import as_vector, as_matrix, diag, inner
from ufl.operators import ListTensor, Zero
from math import sqrt

a: float = (sqrt(3.0)-3.0)/6.0
b: float = (sqrt(3.0)+3.0)/6.0
c: float = -sqrt(3.0)/3.0
d: float = sqrt(2.0)/2.0

E: list[ListTensor] = [diag(as_vector([a,b,c])), diag(as_vector([b,a,c])), as_matrix([[0,d,0],[d,0,0],[0,0,0]]), as_matrix([[0,0,d],[0,0,0],[d,0,0]]), as_matrix([[0,0,0],[0,0,d],[0,d,0]])]

def qvectorfy(tensor: ListTensor | Zero) -> ListTensor | Zero:
    """ qtensor -> qvector """
    return as_vector([inner(E_i,tensor) for E_i in E])

def qtensorfy(vector: ListTensor | Zero) -> ListTensor | Zero:
    """ qvector -> qtensor """
    return sum([vector_i*E_i for vector_i, E_i in zip(vector,E)])
