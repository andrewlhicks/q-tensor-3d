""" The purpose of this module is to add additional functionality to sympy.
The functions here are intended to be used on sympy objects only. """

from sympy import *

# Basis for the 3D Q-tensor

a = (sqrt(3.0)-3.0)/6.0
b = (sqrt(3.0)+3.0)/6.0
c = -sqrt(3.0)/3.0
d = sqrt(2.0)/2.0

E = [diag(a,b,c), diag(b,a,c), Matrix([[0,d,0],[d,0,0],[0,0,0]]), Matrix([[0,0,d],[0,0,0],[d,0,0]]), Matrix([[0,0,0],[0,0,d],[0,d,0]])]

# Test

class DummyMat(MutableDenseMatrix):
    def __new__(cls): 
        return super(QTensor,cls).__new__(cls,zeros(3,3))
    def dx(self,dim_no):
        return tensorfy(self.vect.dx(dim_no))

def arrange(Q,permutation):
    temp_mat = [zeros(3,3),zeros(3,3),zeros(3,3)]
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                temp_mat[kk][ii,jj] += Q.dx(kk)[permutation[ii],permutation[jj]]
    return temp_mat

def einstein(Q,P,permQ,permP):
    Q_arr = arrange(Q,permQ)
    P_arr = arrange(P,permP)
    summ = 0
    for ii in range(3):
        summ += innerp(Q_arr[ii],P_arr[ii])
    return summ

def checkIfLinear(expression,*args):
    for arg in args:
        if not isinstance(arg,QVector) and not isinstance(arg,AbstractVectorGradient):
            raise TypeError('Must be checking for linearity with respect to type QVector or AbstractVectorGradient.')
    for arg in args:
        if isinstance(arg,QVector):
            for ii in range(5):
                secondderivative = diff(expression,arg[ii],2)
                if not secondderivative.is_zero:
                    return False
            for ii in range(5):
                expression = expression.subs(arg[ii],0)
            if not expression.is_zero:
                return False
        elif isinstance(arg,AbstractVectorGradient):
            for ii in range(5):
                for jj in range(3):
                    secondderivative = diff(expression,arg[ii,jj],2)
                    if not secondderivative.is_zero:
                        return False
            for ii in range(5):
                for jj in range(3):
                    expression = expression.subs(arg[ii,jj],0)
            if not expression.is_zero:
                return False
    return True

class BilinearForm:
    def __init__(self,*args):
        assert(len(args)>=4)
        if not isinstance(args[len(args)-4],AbstractVectorGradient):
            raise TypeError('Penultimate arg must be type AbstractVectorGradient.')
        if not isinstance(args[len(args)-3],QVector):
            raise TypeError('Final arg must be type QVector.')
        if not isinstance(args[len(args)-2],AbstractVectorGradient):
            raise TypeError('Penultimate arg must be type AbstractVectorGradient.')
        if not isinstance(args[len(args)-1],QVector):
            raise TypeError('Final arg must be type QVector.')
        
        self.Dtrialfunc = args[len(args)-4]
        self.trialfunc = args[len(args)-3]
        self.Dtestfunc = args[len(args)-2]
        self.testfunc = args[len(args)-1]

        self.forms = []

        for arg in args:
            self.forms.append(arg)
    
        for ii in range(0,4):
            del self.forms[-1]

        for form in self.forms:
            if (not checkIfLinear(self.trialfunc) and not checkIfLinear(self.Dtrialfunc)) or (not checkIfLinear(self.trialfunc) and not checkIfLinear(self.Dtrialfunc)):
                raise ValueError('All forms must be linear.')
    def __call__(self):
        expr = 0
        for ii in range(0,len(self.forms)):
            expr += self.forms[ii]
        return uflfy(expr)
    def add_form(self,*addenda):
        for addendum in addenda:
            if not checkIfLinear(addendum,self.trialfunc,self.testfunc) and not checkIfLinear(addendum,self.Dtrialfunc,self.Dtestfunc):
                raise ValueError('All forms added must be linear.')
            self.forms.append(addendum)

class LinearForm:
    def __init__(self,*args):
        assert(len(args)>=2)
        if not isinstance(args[len(args)-2],AbstractVectorGradient):
            raise TypeError('Penultimate arg must be type AbstractVectorGradient.')
        if not isinstance(args[len(args)-1],QVector):
            raise TypeError('Final arg must be type QVector.')
        
        self.Dtestfunc = args[len(args)-2]
        self.testfunc = args[len(args)-1]

        self.forms = []

        for arg in args:
            self.forms.append(arg)
    
        for ii in range(0,2):
            del self.forms[-1]

        for form in self.forms:
            if not checkIfLinear(form,self.Dtestfunc) and not checkIfLinear(form,self.testfunc):
                raise ValueError('All forms must be linear.')
    def __call__(self):
        expr = 0
        for ii in range(0,len(self.forms)):
            expr += self.forms[ii]
        return uflfy(expr)
    def add_form(self,*addenda):
        for addendum in addenda:
            if not checkIfLinear(addendum,self.Dtestfunc) and not checkIfLinear(addendum,self.testfunc):
                raise ValueError('All forms added must be linear.')
            self.forms.append(addendum)

# Operations

def innerp(A,B):
    """ Returns the Frobenius inner product of two matrices. """
    return trace(A.T*B)

def outerp(A,B):
    """ Returns the "outer" or "tensor" product of two matrices. """
    return A*B.T

# Other

def checkIfQTensor(obj):
    """ Checks if 'obj' is a Q-tensor. I want to be able to check for
    tracelessness and symmetry, but in practice sometimes errors come up when
    attempting this. """

    if not isinstance(obj,MutableDenseMatrix):
        raise TypeError('Must be type MutableDenseMatrix.')
    elif not obj.shape == (3,3):
        raise ShapeError('Shape must be (3, 3).')
    elif not obj.is_symmetric():
        raise ValueError('Must be symmetric.')
    # The following code does not work as intended:
    # elif not trace(obj).is_zero:
    #     raise ValueError('Must have trace 0')

def checkIfVector(obj,dim):
    """ Checks if 'obj' is a Sympy vector of dimension 'dim'. """
    if not isinstance(obj,MutableDenseMatrix):
        raise TypeError('Must be type MutableDenseMatrix.')
    elif not obj.shape == (dim,1):
        raise ShapeError(f'Shape must be ({dim}, 1).')

def tensorfy(vector):
    """ Returns the Q-Tensor form of any 5D vector. """

    checkIfVector(vector,5)

    tensor = zeros(3,3)

    for ii in range(5):
        tensor += vector[ii]*E[ii]

    return tensor

def uflfy(expression):
    """ Returns the UFL code for a scalar or a QVector. First checks if
    'expression' is a matrix. If not, then returns C code for the expression. This
    is a crude way to check for a scalar, but in practice it should work.
    Otherwise, checks if the expression is a 5D vector and returns the C code for
    that expression. """

    if not isinstance(expression,MutableDenseMatrix):
        return ccode(expression)
    else:
        checkIfVector(expression,5)
        return 'as_vector([' + ','.join([ccode(expression[ii]) for ii in range(5)]) + '])'

def variationalDerivative(lagrangian,*args):
    """ Given an instance of Lagrangian, returns a GeneralForm of order 2 which
    represents the first variational derivative of the Lagrangian object respect
    to the last four objects, i.e. two parameters and their derivatives. """

    if not isinstance(lagrangian,Lagrangian):
        raise TypeError('First positional argument must be type Lagrangian.')

    args = list(args)

    if len(args) != 4:
        raise TypeError('Must have exactly 2 arguments per parameter, that is, 4 arguments.')

    for ii in range(0,4,2):
        if not isinstance(args[ii],AbstractVectorGradient):
            raise TypeError('Odd numbered arguments must be type AbstractVectorGradient.')
        if not isinstance(args[ii+1],QVector):
            raise TypeError('Even numbered arguments must be type QVector.')

    ##########################################################################################################    

    # Here, 'ph' stands for "placeholder"

    Dph = [args[ii] for ii in range(0,4,2)]
    ph = [args[ii] for ii in range(1,4,2)]

    tau = Symbol('tau')
    expr = lagrangian(Dph[0]+tau*Dph[1],ph[0]+tau*ph[1])
    expr = diff(expr,tau).subs(tau,0)
    
    return GeneralForm(expr,Dph[0],ph[0],Dph[1],ph[1])

def secondVariationalDerivative(binaryform,*args):
    """ Given an instance of a GeneralForm of order 2, returns the variational
    derivative as a GeneralForm of order 3. """

    if not isinstance(binaryform,GeneralForm):
        raise TypeError('First positional argument must be type GeneralForm.')
    elif not binaryform.order == 2:
        raise ValueError('GeneralForm must be of order 2.')

    args = list(args)
    
    if len(args) != 6:
        raise TypeError('Must have exactly 2 arguments per parameter, that is, 6 arguments.')

    for ii in range(0,6,2):
        if not isinstance(args[ii],AbstractVectorGradient):
            raise TypeError('Odd numbered arguments must be type AbstractVectorGradient.')
        if not isinstance(args[ii+1],QVector):
            raise TypeError('Even numbered arguments must be type QVector.')

    ##########################################################################################################

    # Here, 'ph' stands for "placeholder"

    Dph = [args[ii] for ii in range(0,len(args),2)]
    ph = [args[ii] for ii in range(1,len(args),2)]

    tau = Symbol('tau')
    expr = binaryform(Dph[0]+tau*Dph[1],ph[0]+tau*ph[1],Dph[2],ph[2])
    expr = diff(expr,tau).subs(tau,0)

    return GeneralForm(expr,Dph[0],ph[0],Dph[1],ph[1],Dph[2],ph[2])

def vectorfy(tensor): 
    """ Returns the vector form of a Q-tensor. Checks if 'tensor' is a Q-tensor in
    the mathematical sense, then returns the corresponding vector. """

    checkIfQTensor(tensor)

    vector = zeros(5,1)

    for ii in range(5):
        vector[ii] += innerp(tensor,E[ii])

    return vector

# Classes

class UserDefinedFunction:
    """ Defines a function that can be inputted by the user in a YAML settings
    file simply by entering a string in Sympy format. The spatial variables must
    be 'x0', 'x1', and 'x2'. """
    def __init__(self,func_string):
        self.func_string = func_string
    def  __call__(self):
        x0, x1, x2 = symbols('x0 x1 x2')
        return eval(self.func_string)

class AbstractVectorGradient(Matrix):
    """ Defines a gradient or Jacobian matrix for a given AbstractVector, using
    the built in dx() method of the AbstractVector. """
    def __new__(cls,abstractvector,dim=3):
        if not isinstance(abstractvector,AbstractVector):
            raise TypeError('Must be type AbstractVector.')

        abstractvectorgradient = Matrix([])

        for ii in range(dim):
            abstractvectorgradient = abstractvectorgradient.col_insert(ii,abstractvector.dx(ii))
        
        return super(AbstractVectorGradient,cls).__new__(cls,abstractvectorgradient)

class AbstractVector(Matrix):
    """ Defines a 'dim'-dimensional vector, whose entries are symbols labeled by
    'name' and subscripted from 0 to dim-1. The name chosen should match the
    variable name that will later be used by Firedrake. The dx() method is a
    similar vector but with the '.dx()' suffix attached to each entry of the
    vector, to indicate in UFL code that the derivative is being taken.
    """

    def __new__(cls,name,dim=3):
        if not isinstance(name,str):
            raise TypeError('Must be a string.')

        vector = zeros(dim,1)

        for ii in range(dim):
            vector[ii] += Symbol(f'{name}[{ii}]')
        
        return super(AbstractVector,cls).__new__(cls,vector)

    def __init__(self,name,dim=3):
        self.name = name
        self.dim = dim
        
        self.grad = AbstractVectorGradient(self)
    
    def dx(self,dim_no):
        vector = zeros(self.dim,1)

        for ii in range(self.dim):
            vector[ii] += Symbol(f'{self.name}[{ii}].dx({dim_no})')

        return vector # Ideally, another AbstractVector would be returned, but in practice this is hard

class Lagrangian:
    """ Defines a Lagrangian, i.e., the integrand of an energy functional, in
    terms of user-inputted 'Dph' and 'ph', i.e. the parameter and its derivative
    that the Langrangian will defined with respect to. """

    def __init__(self,expr,Dph,ph):
        self.expr = expr
        self.Dph = Dph
        self.ph = ph
    def __call__(self,Dq=None,q=None):
        if Dq is None and q is None:
            return self.expr
        elif Dq is None:
            raise TypeError('Lagrangian function missing 1 required positional argument: \'Dq\'.')
        elif q is None:
            raise TypeError('Lagrangian function missing 1 required positional argument: \'q\'.')
        elif not isinstance(Dq,AbstractVectorGradient):
            raise TypeError('Argument \'Dq\' for Lagrangian function must be type AbstractVectorGradient.')
        elif not isinstance(q,QVector):
            raise TypeError('Argument \'q\' for Lagrangian function must be type QVector.')
        
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

class GeneralForm:
    """ Returns a more general form of the Lagrangian. More parameters and their
    derivatives are allowed besides just one. The order of the GeneralForm is the
    number of pairs of parameters and their derivative. """
    
    def __init__(self,expr,*args):
        self.expr = expr
        
        args = list(args)
        
        # If no args are given
        
        if len(args) == 0:
            raise TypeError('Expression must be in terms of parameters; none entered.')
        
        # If args are given
        
        # Error handler
        
        if len(args) % 2 != 0:
            raise TypeError('Must have exactly 2 arguments per parameter, that is, an even number of arguments.')
        self.order = int(len(args)/2)
        for ii in range(0,len(args),2):
            if not isinstance(args[ii],AbstractVectorGradient):
                raise TypeError('Odd numbered arguments must be type AbstractVectorGradient.')
            if not isinstance(args[ii+1],QVector):
                raise TypeError('Even numbered arguments must be type QVector.')
        
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
        
        if len(args) != self.order*2:
            raise TypeError(f'Must have exactly {self.order} parameters with 2 arguments each, that is, {self.order*2} arguments.')
        
        # Function
        
        Dq = [args[ii] for ii in range(0,len(args),2)]
        q = [args[ii] for ii in range(1,len(args),2)]
        

        for ii in range(self.order):
            lagrangian_in_one_var = Lagrangian(expr,self.Dph[ii],self.ph[ii])
            expr = lagrangian_in_one_var(Dq[ii],q[ii])

        return expr

class QTensor(Matrix):
    """ Defines a QTensor given a QVector object. Assigns the QVector to .vect.
    The dx() method is the tensorfied dx() of the QVector. """
    def __new__(cls,qvector):
        if not isinstance(qvector,QVector):
            raise TypeError('Must be type QVector.')
        return super(QTensor,cls).__new__(cls,tensorfy(qvector))
    def __init__(self,qvector):
        self.vect = qvector
    def dx(self,dim_no):
        return tensorfy(self.vect.dx(dim_no))

class QVector(AbstractVector):
    """ Defines an AbstractVector of dimension 5. Adds a .tens variable which is
    the QTensor corresponding to the QVector. """
    def __new__(cls,name):
        return super(QVector,cls).__new__(cls,name,5)
    def __init__(self,name):
        super().__init__(name,5)
        self.tens = QTensor(self)

# END OF CODE