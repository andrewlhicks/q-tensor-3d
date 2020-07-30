from sympy import sqrt, diag, Matrix

# Basis for the 3D Q-tensor

a = (sqrt(3.0)-3.0)/6.0
b = (sqrt(3.0)+3.0)/6.0
c = -sqrt(3.0)/3.0
d = sqrt(2.0)/2.0

E = [diag(a,b,c), diag(b,a,c), Matrix([[0,d,0],[d,0,0],[0,0,0]]), Matrix([[0,0,d],[0,0,0],[d,0,0]]), Matrix([[0,0,0],[0,0,d],[0,d,0]])]

# Operations

def innerp(A,B): # Returns the Frobenius inner product of two matrices
    from sympy import trace
    return trace(A.T*B)

def outerp(A,B): # Returns the outer or tensor product of two matrices
    return A*B.T

# Other

def isMatrix(obj): # Checks to see if the object is a MutableDenseMatrix. If so, the returns the Matrix shape. If not, returns 0.
    from sympy import MutableDenseMatrix

    if isinstance(obj,MutableDenseMatrix):
        return obj.shape
    else:
        return 0

def uflfy(expression): # Checks to see if the sympy expression is a scalar or a 5D vector, then returns the UFL code in string form
    from sympy import ccode

    # Test if expression is scalar

    if isMatrix(expression) == 0:
        return ccode(expression)
    
    # Test if expression is vector

    elif isMatrix(expression)[1] == 1:
        # Check if 5D

        if isMatrix(expression)[0] == 5:
            return f"as_vector([{ccode(expression[0])},{ccode(expression[1])},{ccode(expression[2])},{ccode(expression[3])},{ccode(expression[4])}])"

        # Else raise a ValueError

        else:
            raise ValueError("Vector must be 5D.")

    # Else raise a ValueError

    else:
        raise ValueError("Must be scalar or vector.")

def vectorfy(tensor): # Returns the vector form of a Q-tensor
    if isMatrix(tensor) == (3,3):
        return Matrix([innerp(tensor,E[0]),innerp(tensor,E[1]),innerp(tensor,E[2]),innerp(tensor,E[3]),innerp(tensor,E[4])])
    else:
        raise ValueError('Must be a 3x3 tensor')

# Classes

class UserDefinedFunction:
    def __init__(self,func_string):
        self.func_string = func_string
    def  __call__(self):
        import sympy as sp
        x0, x1, x2 = sp.symbols('x0 x1 x2')
        return eval(self.func_string)

class Vector(Matrix):
    def __new__(cls,name,dim=3):
        from sympy import symbols
        vector = []
        for ii in range(dim):
            vector.append(symbols(f'{name}[{ii}]'))
        
        return super(Vector,cls).__new__(cls,vector)
    def __init__(self,name,dim=3):
        self.name = name # For 'name', choose the variable name that Firedrake will later use
        self.dim = 5

class QVector(Vector):
    def __new__(cls,name):
        return super(QVector,cls).__new__(cls,name,5)
    def __init__(self,name):
        super().__init__(name,5)

        from sympy import symbols, zeros

        grad = []
        grad_row = []

        for ii in range(5):
            for jj in range(3):
                grad_row.append(symbols(f'{name}[{ii}].dx({jj})'))
            grad.append(grad_row)
            grad_row = []

        self.grad = Matrix(grad)

        tensor = zeros(3,3)

        for ii in range(5):
            tensor += self[ii]*E[ii]

        self.ten = tensor

# END OF CODE