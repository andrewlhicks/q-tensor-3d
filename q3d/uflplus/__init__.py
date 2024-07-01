from ufl import *
from ufl.classes import * # allows for evaluating UFL cache strings into UFL objects
from firedrake.ufl_expr import * # overwrites some of ufl.classes using Firedrake-specific UFL object defintinions. Needed for Firedrake solvers to work properly
from ufl.operators import ListTensor, Zero # allows for type hinting/checking for UFL objects

def vectorfy(tensor: ListTensor | Zero) -> ListTensor | Zero:
    """ qtensor -> qvector """
    from math import sqrt

    a: float = (sqrt(3.0)-3.0)/6.0
    b: float = (sqrt(3.0)+3.0)/6.0
    c: float = -sqrt(3.0)/3.0
    d: float = sqrt(2.0)/2.0

    E: list[ListTensor] = [diag(as_vector([a,b,c])), diag(as_vector([b,a,c])), as_matrix([[0,d,0],[d,0,0],[0,0,0]]), as_matrix([[0,0,d],[0,0,0],[d,0,0]]), as_matrix([[0,0,0],[0,0,d],[0,d,0]])]

    m: ListTensor = as_vector([inner(E_i,tensor) for E_i in E])
    return m

def smooth_transition(x, *, I: list[int | float]):
    """ A smooth transition function in the variable x on the
    interval I = [a,b], with the value of 0 at x=a and 1 at x=b """

    a = I[0]
    b = I[1]

    # smooth function f defined by the exponential
    def f(x):
        return conditional(gt(x,0), exp(-1/x), 0)
    # smooth function g such that g(x)=0 for x<=0 and g(x)=1 for x>=1
    def g(x):
        return f(x)/(f(x)+f(1-x))
    
    return g((x-a)/(b-a))