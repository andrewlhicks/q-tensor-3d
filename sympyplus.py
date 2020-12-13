""" The purpose of this module is to add additional functionality to sympy.
The functions here are intended to be used on sympy objects only. """

from sympy import *

# Set up variable with respect to which we will take derivatives

x = [Symbol('x0'),Symbol('x1'),Symbol('x2')]

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

def levi_civita(i,j,k):
    indices = [i,j,k]

    for index in indices:
        if not (index == 0 or index == 1 or index == 2):
            raise TypeError(f'Index must be 0, 1, or 2; {index} was given.')

    def index_rearrange(indices):
        if indices[0] == indices[1] or indices[1] == indices[2] or indices[0] == indices[2]:
            return [0,0,0]
        elif indices[0] == 0:
            # print(f'Returing indices: {indices}')
            return indices
        else:
            new_indices = [0,0,0]
            for i in range(3):
                new_indices[i-1] += indices[i]
            # print(f'Old indices: {indices}, New indices: {new_indices}')
            return index_rearrange(new_indices)

    indices = index_rearrange(indices)
    # print(indices)

    if indices == [0,0,0]:
        return 0
    elif indices == [0,1,2]:
        return 1
    elif indices == [0,2,1]:
        return -1
    else:
        raise ValueError('Indices values messed up.')

def mixedp(A,B):
        """ Returns the mixed product of QTensor B with the derivative of QTensor A, that is, epsilon(i,j,k)A(l,j,k)B(i,j) """
        product = 0

        for ii in range(3):
            for jj in range(3):
                for kk in range(3):
                    for ll in range(3):
                        product += levi_civita(ii,kk,ll)*A.dx(kk)[ll,jj]*B[ii,jj]

        return product

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

def isLinear(expression,*args):
    for arg in args:
        if not isinstance(arg,QVector) and not isinstance(arg,AbstractVectorGradient):
            raise TypeError('Must be checking for linearity with respect to type QVector or AbstractVectorGradient.')
    for arg in args:
        if isinstance(arg,QVector):
            for ii in range(5):
                secondderivative = diff(expression,arg[ii],2)
                if not secondderivative.is_zero:
                    # print('Second der not zero')
                    return False
            for ii in range(5):
                expression = expression.subs(arg[ii],0)
            if not expression.is_zero:
                # print('0 maps to nonzero')
                return False
        elif isinstance(arg,AbstractVectorGradient):
            for ii in range(5):
                for jj in range(3):
                    secondderivative = diff(expression,arg[ii,jj],2)
                    if not secondderivative.is_zero:
                        # print('Second der not zero')
                        return False
            for ii in range(5):
                for jj in range(3):
                    expression = expression.subs(arg[ii,jj],0)
            if not expression.is_zero:
                # print('0 maps to nonzero')
                return False
    return True

class Param:
    def __init__(self,param):
        if isinstance(param,Param):
            self.der = param.der
            self.vec = param.vec
        else:
            if not isinstance(param,list):
                raise TypeError('Parameters must be lists.')
            elif len(param) != 2:
                raise ValueError('Parameter must be be a list of length 2.')
            elif not isinstance(param[0],AbstractVectorGradient):
                raise TypeError('First argument of parameter must be type AbstractVectorGradient.')
            elif not isinstance(param[1],QVector):
                raise TypeError('Second argument of parameter must be type QVector.')
            self.der = param[0]
            self.vec = param[1]

def isLinearParam(expression,param):
    if not isinstance(param,Param):
        raise TypeError('Second positional argument must be type Param.')

    if isLinear(expression,param.der) or isLinear(expression,param.vec):
        return True
    else:
        return False

def checkIfParam(param):
    if not isinstance(param,list):
        raise TypeError('Parameters must be lists.')
    elif len(param) != 2:
        raise ValueError('Parameter must be be a list of length 2.')
    elif not isinstance(param[0],AbstractVectorGradient):
        raise TypeError('First argument of parameter must be type AbstractVectorGradient.')
    elif not isinstance(param[1],QVector):
        raise TypeError('Second argument of parameter must be type QVector.')

class lhsForm:
    def __init__(self,trial_func,test_func,name=None,forms=[]):
        if not isinstance(forms,list):
            raise TypeError('\'forms\' must be a List of items of type GeneralForm.')
        self.trial_func = Param(trial_func)
        self.test_func = Param(test_func)
        self.name = name
        self.forms = []
        self.add_form(*forms)

    def __call__(self):
        expr = 0
        for form in self.forms:
            if not isLinearParam(form.expr,self.trial_func):
                # raise ValueError(f'The form \'{form}\' of \'{self}\' is nonlinear in the trial function.')
                print(f'Sympyplus detects form \'{form}\' of \'{self}\' to be nonlinear in the trial function; however, this is an error.')
            elif not isLinearParam(form.expr,self.test_func):
                # raise ValueError(f'The form \'{form}\' of \'{self}\' is nonlinear in the test function.')
                print(f'Sympyplus detects form \'{form}\' of \'{self}\' to be nonlinear in the test function; however, this is an error.')
            expr += form.expr
        return uflfy(expr)

    def __repr__(self):
        return 'Untitled lhsForm' if self.name == None else f'lhsForm {self.name}'

    def add_form(self,*forms):
        for form in forms:
            if not isinstance(form,GeneralForm):
                raise TypeError('Form must be type GeneralForm.')
            self.forms.append(form)

class rhsForm:
    def __init__(self,test_func,name=None,forms=[]):
        if not isinstance(forms,list):
            raise TypeError('\'forms\' must be a List of items of type GeneralForm.')
        self.test_func = Param(test_func)
        self.name = name
        self.forms = []
        self.add_form(*forms)

    def __call__(self):
        expr = 0
        for form in self.forms:
            if not isLinearParam(form.expr,self.test_func):
                # raise ValueError(f'The form \'{form}\' of \'{self}\' is nonlinear.')
                print(f'Sympyplus detects form \'{form}\' of \'{self}\' to be nonlinear; however, this is an error.')
            expr += form.expr
        return uflfy(expr)

    def __repr__(self):
        return 'Untitled rhsForm' if self.name == None else f'rhsForm {self.name}'

    def add_form(self,*forms):
        for form in forms:
            if not isinstance(form,GeneralForm):
                raise TypeError('Form must be type GeneralForm.')
            self.forms.append(form)

def newtonsMethod(lhs_form,rhs_form,const_func,trial_func,test_func):
    const_func = Param(const_func)
    trial_func = Param(trial_func)
    test_func = Param(test_func)

    new_lhs_form = lhsForm(trial_func,test_func)
    for form in lhs_form.forms:
        new_lhs_form.add_form(secondVariationalDerivative(form,const_func,trial_func,test_func))

    new_rhs_form = rhsForm(test_func)
    for form in lhs_form.forms:
        new_rhs_form.add_form(form.mul(-1).new_params(const_func,test_func))
    for form in rhs_form.forms:
        new_rhs_form.add_form(form.new_params(test_func))

    return (new_lhs_form,new_rhs_form)


# Operations

def innerp(A,B):
    """ Returns the Frobenius inner product of two matrices. """
    return trace(A.T*B)

def outerp(A,B):
    """ Returns the "outer" or "tensor" product of two matrices. """
    return A*B.T

def fnorm(A):
    """ Returns the Frobenius norm of the matrix. """
    return sqrt(innerp(A,A))

# Other

def checkIfQTensor(obj):
    """ Checks if 'obj' is a Q-tensor. I want to be able to check for
    tracelessness and symmetry, but in practice sometimes errors come up when
    attempting this. """

    if not isinstance(obj,MutableDenseMatrix):
        raise TypeError('Must be type MutableDenseMatrix.')
    elif not obj.shape == (3,3):
        raise ShapeError('Shape must be (3, 3).')
    # elif not obj.is_symmetric():
    #     raise ValueError('Must be symmetric.')
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

def variationalDerivative(lagrangian,*params,name=None):
    """ Given an instance of Lagrangian, returns a GeneralForm of order 2 which
    represents the first variational derivative of the Lagrangian object respect
    to the last four objects, i.e. two parameters and their derivatives. """

    if not isinstance(lagrangian,GeneralForm):
        raise TypeError('First positional argument must be type GeneralForm.')
    elif lagrangian.order != 1:
        raise ValueError('GeneralForm must be of order 1.')

    if len(params) != 2:
        raise TypeError('Must have exactly 2 parameters.')

    params = [Param(param) for param in params]

    ##########################################################################################################    

    # Compute the derivative

    tau = Symbol('tau')
    expr = lagrangian([params[0].der+tau*params[1].der,params[0].vec+tau*params[1].vec])
    expr = diff(expr,tau).subs(tau,0)
    
    return GeneralForm(expr,[params[0].der,params[0].vec],[params[1].der,params[1].vec],name=name)

def secondVariationalDerivative(binaryform,*params,name=None):
    """ Given an instance of a GeneralForm of order 2, returns the variational
    derivative as a GeneralForm of order 3. """

    if not isinstance(binaryform,GeneralForm):
        raise TypeError('First positional argument must be type GeneralForm.')
    elif not binaryform.order == 2:
        raise ValueError('GeneralForm must be of order 2.')
    
    if len(params) != 3:
        raise TypeError('Must have exactly 3 parameters.')
    
    params = [Param(param) for param in params]

    ##########################################################################################################

    tau = Symbol('tau')
    expr = binaryform([params[0].der+tau*params[1].der,params[0].vec+tau*params[1].vec],[params[2].der,params[2].vec])
    expr = diff(expr,tau).subs(tau,0)

    if name == None:
        name = f'Der of {binaryform.name}'

    return GeneralForm(expr,[params[0].der,params[0].vec],[params[1].der,params[1].vec],[params[2].der,params[2].vec],name=name)

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

class GeneralForm:
    """ Returns a more general form of the Lagrangian. More parameters and their
    derivatives are allowed besides just one. The order of the GeneralForm is the
    number of pairs of parameters and their derivative. """
    
    def __init__(self,expr,*params,name=None):
        self.expr = expr
        self.name = name
        self.order = len(params)

        if self.order == 0:
            raise TypeError('Expression must be in terms of parameters; none entered.')
        
        self.params = [Param(param) for param in params]

    def __call__(self,*new_params):
        if len(new_params) == 0:
            return self.expr
        elif len(new_params) != self.order:
            raise ValueError(f'The number of params must be equal to the order of the form ({self.order}). {len(new_params)} were given.')
        else:
            return self.new_params(*new_params).expr

    def __repr__(self):
        if self.name is not None:
            return self.name
        else:
            return 'Unnamed'

    def rename(self,name):
        if not isinstance(name,str):
            raise TypeError('Name must be a string.')
        self.name = name

    def new_params(self,*new_params):
        if len(new_params) != self.order:
            raise ValueError(f'The number of new params must be equal to the order of the form ({self.order}). {len(new_params)} were given.')
        
        new_params = [Param(new_param) for new_param in new_params]

        # Initialize

        new_expr = self.expr

        # Function

        for kk in range(self.order):            
            old_param = self.params[kk]
            new_param = new_params[kk]

            # Substitute old parameter derivative for new parameter derivative

            for jj in range(3):
                for ii in range(5):
                    old = old_param.der[ii,jj]
                    new = new_param.der[ii,jj]
                    new_expr = new_expr.subs(old,new)
            
            # Substitute old parameter vector for new parameter vector

            for ii in range(5):
                old = old_param.vec[ii]
                new = new_param.vec[ii]
                new_expr = new_expr.subs(old,new)

        return GeneralForm(new_expr,*new_params,name=self.name)

    def mul(self,n):
        new_expr = self.expr * n
        return GeneralForm(new_expr,*self.params,name=self.name)

    def checkParam(self,param):
        if not isinstance(param,list):
                raise TypeError('Parameters must be lists.')
        elif len(param) != 2:
            raise ValueError('Parameter must be be a list of length 2.')
        elif not isinstance(param[0],AbstractVectorGradient):
            raise TypeError('First argument of parameter must be type AbstractVectorGradient.')
        elif not isinstance(param[1],QVector):
            raise TypeError('Second argument of parameter must be type QVector.')

    def uflfy(self):
        return uflfy(self.expr)

class Lagrangian(GeneralForm):
    def __init__(self,expr,*params,name=None):
        if len(params) != 1:
            raise ValueError(f'Lagrangian must contain only one parameter. {len(params)} were given.')
        return super.__init__(self,expr,*params,name=None)

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

# Set up Qvector objects

nu = AbstractVector('nu')

q = QVector('q')
Dq = q.grad
Q = q.tens

p = QVector('p')
Dp = p.grad
P = p.tens

qp = QVector('q_prev')
Dqp = qp.grad
QP = qp.tens

qnp = QVector('q_newt_prev')
Dqnp = qnp.grad
QNP = qnp.tens

f = QVector('f')
f_gam = QVector('f_gam')

# END OF CODE