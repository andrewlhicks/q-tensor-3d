from firedrake import *
import sympy as sp

a = (sp.sqrt(3.0)-3.0)/6.0
b = (sp.sqrt(3.0)+3.0)/6.0
c = -sp.sqrt(3.0)/3.0
d = sp.sqrt(2.0)/2.0

E = [sp.diag(a,b,c), sp.diag(b,a,c), sp.Matrix([[0,d,0],[d,0,0],[0,0,0]]), sp.Matrix([[0,0,d],[0,0,0],[d,0,0]]), sp.Matrix([[0,0,0],[0,0,d],[0,d,0]])]

def isMatrix(obj):
    if isinstance(obj,sp.MutableDenseMatrix):
        return obj.shape
    else:
        return 0

def uflfy(expression):
    if isMatrix(expression) == 0:                           # Test if expression is scalar
        return sp.ccode(expression)
    elif isMatrix(expression)[1] == 1:                      # Test if expression is vector
        if isMatrix(expression)[0] == 2:                        # Check if 2D
            return f"as_vector([{sp.ccode(expression[0])},{sp.ccode(expression[1])}])"
        elif isMatrix(expression)[0] == 5:                      # Check if 3D
            return f"as_vector([{sp.ccode(expression[0])},{sp.ccode(expression[1])},{sp.ccode(expression[2])},{sp.ccode(expression[3])},{sp.ccode(expression[4])}])"
        else:
            raise ValueError("Vector must be 2D or 5D.")
    else:
        raise ValueError("Must be scalar or vector.")

def innerp(A,B):
    return sp.trace(A.T*B)

def outerp(A,B):
    return A*B.T

def vectorfy(tensor):
    if isMatrix(tensor) == (2,2):
        return sp.Matrix([innerp(tensor,E[0]),innerp(tensor,E[1])])
    if isMatrix(tensor) == (3,3):
        return sp.Matrix([innerp(tensor,E[0]),innerp(tensor,E[1]),innerp(tensor,E[2]),innerp(tensor,E[3]),innerp(tensor,E[4])])
    else:
        raise ValueError('Must be a 2x2 or 3x3 tensor')