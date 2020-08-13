from sympy import *

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

def variationalDerivative(energy,Dparam1,param1,Dparam2,param2):
    tau = Symbol('tau')
    expr = energy(param1.grad+tau*param2.grad,param1+tau*param2)
    expr = diff(expr,tau)
    expr = expr.subs(tau,0)
    return BinaryForm(expr,Dparam1,param1,Dparam2,param2)

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

class AbstractVectorGradient(Matrix):
    def __new__(cls,abstract_vector,dim=3):
        gradient_vector = Matrix([])
        for ii in range(dim):
            gradient_vector = gradient_vector.col_insert(ii,abstract_vector.dx(ii))
        
        return super(AbstractVectorGradient,cls).__new__(cls,gradient_vector)

class AbstractVector(Matrix):
    def __new__(cls,name,dim=3):
        from sympy import Symbol
        vector = []
        for ii in range(dim):
            vector.append(Symbol(f'{name}[{ii}]'))
        
        return super(AbstractVector,cls).__new__(cls,vector)
    def __init__(self,name,dim=3):
        self.name = name # For 'name', choose the variable name that Firedrake will later use
        self.dim = 5
        
        self.grad = AbstractVectorGradient(self)
    
    def dx(self,dim_no):
        from sympy import Symbol, zeros
        vector = []
        for ii in range(self.dim):
            vector.append(Symbol(f'{self.name}[{ii}].dx({dim_no})'))
        
        return Matrix(vector) # Maybe in the future try to retun another Abstract Vector? But to do this I would need to make AbstractVector more general in terms of name input

class Lagrangian:
    def __init__(self,expr,Dph,ph): # Later on, maybe try to apply this code to the class UserDefinedFunction
        self.expr = expr
        self.Dph = Dph
        self.ph = ph
    def __call__(self,Dq=None,q=None):
        if Dq is None and q is None:
            return self.expr
        elif Dq is None:
            raise TypeError('Lagrangian function missing 1 required positional argument: \'Dq\'')
        elif q is None:
            raise TypeError('Lagrangian function missing 1 required positional argument: \'q\'')
        elif not isinstance(Dq,AbstractVectorGradient):
            raise TypeError('Argument \'Dq\' for Lagrangian function must be type AbstractVectorGradient')
        elif not isinstance(q,QVector):
            raise TypeError('Argument \'q\' for Lagrangian function must be type QVector')
        
        expr = self.expr
        Dph = self.Dph
        ph = self.ph
        
        for jj in range(3):
            for ii in range(5):
                old = Dph[ii,jj]
                new = Dq[ii,jj]
                expr = expr.subs(old,new)
        
        for ii in range(5):
            old = ph[ii]
            new = q[ii]
            expr = expr.subs(old,new)
        
        return expr

class BinaryForm:
    def __init__(self,expr,Dph1,ph1,Dph2,ph2):
        self.expr = expr
        self.Dph1 = Dph1
        self.ph1 = ph1
        self.Dph2 = Dph2
        self.ph2 = ph2
    def __call__(self,Dq=None,q=None,Dp=None,p=None):
        expr = self.expr
        lagrangian_in_var1 = Lagrangian(expr,self.Dph1,self.ph1)
        expr = lagrangian_in_var1(Dq,q)
        lagrangian_in_var2 = Lagrangian(expr,self.Dph2,self.ph2)
        expr = lagrangian_in_var2(Dp,p)
        return expr

class GeneralForm:
    def __init__(self,expr,*args):
        self.expr = expr
        
        args = list(args)
        
        # If no args are given
        
        if len(args) == 0:
            raise TypeError('Expression must be in terms of parameters; none entered')
        
        # If args are given
        
        # Error handler
        
        if len(args) % 2 != 0:
            raise TypeError('Must have exactly 2 arguments per parameter, that is, an even number of arguments')
        self.no_params = int(len(args)/2)
        for ii in range(0,len(args),2):
            if not isinstance(args[ii],AbstractVectorGradient):
                raise TypeError('Odd numbered arguments must be type AbstractVectorGradient')
            if not isinstance(args[ii+1],QVector):
                raise TypeError('Even numbered arguments must be type QVector')
        
        # Function
        
        self.Dph = [args[ii] for ii in range(0,len(args),2)]
        self.ph = [args[ii] for ii in range(1,len(args),2)]
    def __call__(self,*args):
        expr = self.expr
        
        args = list(args)
        
        # If no args are given
        
        if len(args) == 0:
            return expr
        
        # If args are given
        
        # Error handler
        
        if len(args) != self.no_params*2:
            raise TypeError(f'Must have exactly {self.no_params} parameters with 2 arguments each, that is, {self.no_params*2} arguments')
        
        # Function
        
        Dq = [args[ii] for ii in range(0,len(args),2)]
        q = [args[ii] for ii in range(1,len(args),2)]
        

        for ii in range(self.no_params):
            lagrangian_in_one_var = Lagrangian(expr,self.Dph[ii],self.ph[ii])
            expr = lagrangian_in_one_var(Dq[ii],q[ii])

        return expr

class QTensor(Matrix):
    def __new__(cls,q_vector):
        from sympy import zeros
        
        tensor = zeros(3,3)
        for ii in range(5):
            tensor += q_vector[ii]*E[ii]
            
        return super(QTensor,cls).__new__(cls,tensor)
    def __init__(self,q_vector):
        self.vect = q_vector
    def dx(self,dim_no):
        from sympy import zeros

        tensor_dx = zeros(3,3) # perhaps I could make a function to do this?
        for ii in range(5):
            tensor_dx += self.vect.dx(dim_no)[ii]*E[ii]

        return tensor_dx

class QVector(AbstractVector):
    def __new__(cls,name):
        return super(QVector,cls).__new__(cls,name,5)
    def __init__(self,name):
        super().__init__(name,5)

        self.tens = QTensor(self)

# END OF CODE