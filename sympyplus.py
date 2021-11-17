""" The purpose of this module is to add additional functionality to sympy.
The functions here are intended to be used on sympy objects only. """

from typing import Type
from sympy import *

# Basis for the 3D Q-tensor

a = (sqrt(3.0)-3.0)/6.0
b = (sqrt(3.0)+3.0)/6.0
c = -sqrt(3.0)/3.0
d = sqrt(2.0)/2.0

E = [diag(a,b,c), diag(b,a,c), Matrix([[0,d,0],[d,0,0],[0,0,0]]), Matrix([[0,0,d],[0,0,0],[d,0,0]]), Matrix([[0,0,0],[0,0,d],[0,d,0]])]

# FUNCTIONS

# Operations

def innerp(A,B):
    """ Returns the Frobenius inner product of two matrices. """
    return trace(A.T*B)

def outerp(A,B):
    """ Returns the "outer" or "tensor" product of two matrices. """
    return A*B.T

def mixedp(A,B):
    """ Returns the mixed product of QTensor B with the derivative of QTensor A,
    that is, epsilon(i,j,k)A(l,j,k)B(i,j) """

    product = 0

    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                for ll in range(3):
                    product += levi_civita(ii,kk,ll)*A.dx(kk)[ll,jj]*B[ii,jj]

    return product

def levi_civita(i,j,k):
    """ Returns the Levi-Civita symbol evaluated at indices i, j, and k, where
    the indices are in the range [0,1,2]. """

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

# Norms

def fnorm(A):
    """ Returns the Frobenius norm of the matrix. """
    return sqrt(innerp(A,A))

# Checks

def is_linear_param(expression,parameter):
    """ Checks if the expression is linear in this single parameter. """
    if not isinstance(parameter,Param):
        raise TypeError('Parameter given must be type Param.')
    args = parameter.explode()

    # Begin with the second derivative test for all parameters

    for arg in args:
        second_derivative = diff(expression,arg,2)
        if not second_derivative.is_zero:
            return False

    # Next, set all args equal to zero and see if expression becomes zero

    for arg in args:
        expression = expression.subs(arg,0)
    if not expression.is_zero:
        return False

    # If the two tests do not return False, return True

    return True

def checkIfParam(param):
    if not isinstance(param,list):
        raise TypeError('Parameters must be lists.')
    elif len(param) != 2:
        raise ValueError('Parameter must be be a list of length 2.')
    elif not isinstance(param[0],AbstractVectorGradient):
        raise TypeError('First argument of parameter must be type AbstractVectorGradient.')
    elif not isinstance(param[1],QVector):
        raise TypeError('Second argument of parameter must be type QVector.')

def check_istensor(obj):
    """ Raises error if 'obj' is not a Sympy Matrix of shape (3,3). """
    if not isinstance(obj,MutableDenseMatrix):
        raise TypeError('Must be type MutableDenseMatrix.')
    if obj.shape != (3,3):
        raise ShapeError(f'Shape must be (3, 3), not {obj.shape}.')

def istensor(obj):
    """ Checks if 'obj' is a Sympy Matrix of shape (3,3). """
    if not isinstance(obj,MutableDenseMatrix):
        return False
    if obj.shape != (3,3):
        return False
    return True

def check_isvector(obj,dim=None):
    """ Raises error if 'obj' is not a Sympy vector of dimension 'dim'. """
    if not isinstance(obj,MutableDenseMatrix):
        raise TypeError('Must be type MutableDenseMatrix.')
    if obj.shape[1] != 1:
        raise ShapeError(f'Shape must be (*, 1), not {obj.shape}.')
    if dim is not None and obj.shape[0] != dim:
        raise ShapeError(f'Shape must be ({dim}, 1).')

def isvector(obj,dim=None):
    if not isinstance(obj,MutableDenseMatrix):
        return False
    if obj.shape[1] != 1:
        return False
    if dim is not None and obj.shape[0] != dim:
        return False
    return True

# Type alterers

def tensorfy(vector,basis):
    check_isvector(vector)

    if not isinstance(basis,TensorSpaceBasis):
        raise TypeError('Argument "basis" must be type TensorSpaceBasis.')
    if vector.shape[0] != basis.dim:
        raise DimensionError(f'vector and basis must have same dimension, not {vector.shape[0]} and {basis.dim}')

    tensor = zeros(3,3)

    for ii in range(basis.dim):
        tensor += vector[ii]*basis[ii]

    return tensor

def q_tensorfy(vector):
    """ Returns the Q-Tensor form of any 5D vector. """

    check_isvector(vector,5)

    q_basis = TensorSpaceBasis(E)

    return tensorfy(vector,q_basis)

def uflfy(expression):
    """ Returns the UFL code for a scalar or a QVector. First checks if
    'expression' is a matrix. If not, then returns C code for the expression.
    This is a crude way to check for a scalar, but in practice it should work.
    Otherwise, checks if the expression is a 3D or 5D vector and returns the C
    code for that expression. """

    if not isinstance(expression,MutableDenseMatrix):
        return ccode(expression)
    elif isvector(expression,3):
        return 'as_vector([' + ','.join([ccode(expression[ii]) for ii in range(3)]) + '])'
    elif isvector(expression,5):
        return 'as_vector([' + ','.join([ccode(expression[ii]) for ii in range(5)]) + '])'
    else:
        raise TypeError('Must be a vector expression of dimension 3 or 5.')

def vectorfy(tensor):
    """ Returns the vector form of a Q-tensor. Checks if 'tensor' is a Q-tensor
    in the mathematical sense, then returns the corresponding vector. """

    check_istensor(tensor)

    vector = zeros(5,1)

    for ii in range(5):
        vector[ii] += innerp(tensor,E[ii])

    return vector

# Calculus of variations and Newton's method

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
    derivative as a GeneralForm of order 3.

    Note that this in this implementation, we have for a binary form A[Q](P),
    the derivative dA[Q](P,R) with

    Param 0 = Q
    Param 1 = R
    Param 2 = P

    But a better implementation might be

    Param 0 = Q
    Param 1 = P
    Param 2 = R

    This may be fixed later.
    """

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

# CLASSES

class DimensionError(Exception):
    pass

# Vectors and tensors

class AbstractVectorGradient(Matrix):
    """ Defines a gradient or Jacobian matrix for a given AbstractVector, using
    the built in dx() method of the AbstractVector. """

    def __new__(cls,abstractvector,dim=3):
        if not isinstance(abstractvector,AbstractVector):
            raise TypeError('Must be type AbstractVector.')

        abstractvectorgradient = Matrix([])

        for ii in range(dim):
            abstractvectorgradient = abstractvectorgradient.col_insert(ii,abstractvector.dx(ii))

        return super().__new__(cls,abstractvectorgradient)
    
    def __init__(self,abstractvector,dim=3) -> None:
        self.__name = abstractvector.name

    def __repr__(self):
        return f'D{self.name}'
    
    @property
    def name(self):
        return self.__name

class AbstractVector(Matrix):
    """ Defines a 'dim'-dimensional vector, whose entries are symbols labeled by
    'name' and subscripted from 0 to dim-1. The name chosen should match the
    variable name that will later be used by Firedrake. The dx() method is a
    similar vector but with the '.dx()' suffix attached to each entry of the
    vector, to indicate in UFL code that the derivative is being taken. """

    def __new__(cls,name,dim=3):
        if not isinstance(name,str):
            raise TypeError('Must be a string.')

        vector = zeros(dim,1)

        for ii in range(dim):
            vector[ii] += Symbol(f'{name}[{ii}]')

        return super().__new__(cls,vector)

    def __init__(self,name,dim=3):
        self.name = name
        self.dim = dim

        self.grad = AbstractVectorGradient(self)

    def __repr__(self):
        return self.name

    def dx(self,dim_no):
        vector = zeros(self.dim,1)

        for ii in range(self.dim):
            vector[ii] += Symbol(f'{self.name}[{ii}].dx({dim_no})')

        return vector # Ideally, another AbstractVector would be returned, but in practice this is hard

class TensorSpaceBasis(list):
    def __init__(self,basis):
        space_shape = (3,3)

        if not isinstance(basis,list):
            raise TypeError('TensorSpaceBasis must be type "list".')
        for item in basis:
            if not isinstance(item,Matrix):
                raise TypeError('Basis elements must be tye "Matrix".')
            if item.shape != space_shape:
                raise ShapeError(f'Basis elements must have shape {space_shape}')
        self.dim = len(basis)
        list.__init__(self,basis)

class Param:
    """ Defines a Param object from a list of length 2, or a Param object. The
    first item in the list is the first derivative and the second is the
    parameter without derivative. """
    def __init__(self,param):
        if isinstance(param,Param):
            self.der = param.der
            self.vec = param.vec
        else:
            if not isinstance(param,list):
                raise TypeError('Parameters must be lists.')
            elif len(param) != 2:
                raise ValueError('Parameter must be be a list of length 2.')
            self.der = param[0]
            self.vec = param[1]
    def __repr__(self) -> str:
        return f'[{repr(self.der)},{repr(self.vec)}]'
    def __eq__(self,other) -> bool:
        if not isinstance(other,Param):
            raise TypeError('Must compare Param to Param.')
        if self.der == other.der and self.vec == other.vec:
            return True
        else:
            return False
    def explode(self):
        """ Returns the Symbols of the Param as a list. """
        return [self.der[ii,jj] for ii in range(5) for jj in range(3)] + [self.vec[ii] for ii in range(5)]

    # Properties
    @property
    def der(self):
        return self.__der
    @property
    def vec(self):
        return self.__vec
    
    # Setters
    @der.setter
    def der(self,val):
        if not isinstance(val,AbstractVectorGradient):
            raise TypeError('Must be type AbstractVectorGradient.')
        self.__der = val
    @vec.setter
    def vec(self,val):
        if not isinstance(val,QVector):
            raise TypeError('Must be type QVector')
        self.__vec = val

class SymmetricTracelessMatrix(Matrix):
    """ Returns the symmetric traceless part of the matrix. """
    def __new__(cls,array):
        M = Matrix(array)
        M = (1/2)*(M + M.T)
        M = M - trace(M)/3*eye(3)
        return super().__new__(cls,M)

class QTensor(Matrix):
    """ Defines a QTensor given a QVector object. Assigns the QVector to .vect.
    The dx() method is the tensorfied dx() of the QVector. """
    def __new__(cls,qvector):
        if not isinstance(qvector,QVector):
            raise TypeError('Must be type QVector.')
        return super().__new__(cls,q_tensorfy(qvector))
    def __init__(self,qvector):
        self.vect = qvector
    def dx(self,dim_no):
        return q_tensorfy(self.vect.dx(dim_no))

class QVector(AbstractVector):
    """ Defines an AbstractVector of dimension 5. Adds a .tens variable which is
    the QTensor corresponding to the QVector. """
    def __new__(cls,name):
        return super().__new__(cls,name,5)
    def __init__(self,name):
        super().__init__(name,5)
        self.tens = QTensor(self)

# Forms

class GeneralForm:
    """ Returns a more general form of the Lagrangian. More parameters and their
    derivatives are allowed besides just one. The order of the GeneralForm is
    the number of pairs of parameters and their derivative. """

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

class EnergyForm:
    def __init__(self,domain=[],boundary=[]):
        if not isinstance(domain,list) or not isinstance(boundary,list):
            raise TypeError
        for item in domain + boundary:
            if not isinstance(item,GeneralForm):
                raise TypeError
        
        # Initialize params

        params = [[Dq,q],[Dp,p],[Dr,r]] # For now, let the params be fixed, but allow customization later
        self.__params = [Param(param) for param in params]

        # Initialize domain
        self.__domain_0, self.__domain_1, self.__domain_2 = [], [], []
        self.add_domain(*domain)

        # Initialize boundary
        self.__boundary_0, self.__boundary_1, self.__boundary_2 = [], [], []
        self.add_boundary(*boundary)

    def __repr__(self):
        return f'<EnergyForm d={self.domain} b={self.boundary}>'
    
    @property
    def domain(self):
        return [repr(form) for form in self.__domain_0]
    @property
    def boundary(self):
        return [repr(form) for form in self.__boundary_0]

    @property
    def params(self):
        return self.__params

    @property
    def domain_0(self):
        return [form.uflfy() for form in self.__domain_0]
    @property
    def domain_1(self):
        return [form.uflfy() for form in self.__domain_1]
    @property
    def domain_2(self):
        return [form.uflfy() for form in self.__domain_2]
    @property
    def boundary_0(self):
        return [form.uflfy() for form in self.__boundary_0]
    @property
    def boundary_1(self):
        return [form.uflfy() for form in self.__boundary_1]
    @property
    def boundary_2(self):
        return [form.uflfy() for form in self.__boundary_2]

    @property
    def ufl_dict(self):
        return [
            {'domain':self.domain_0,'boundary':self.boundary_0},
            {'domain':self.domain_1,'boundary':self.boundary_1},
            {'domain':self.domain_2,'boundary':self.boundary_2}
        ]
    
    def differentiate_forms(self,*forms):
        # Preliminary checks
        for form in forms:
            if not isinstance(form,GeneralForm):
                raise TypeError('"forms" must be a list of items of type GeneralForm.')
            if form.order != 1:
                raise ValueError
        # Get derivatives of forms
        forms_0 = list(forms)
        forms_1 = [variationalDerivative(form,self.params[0],self.params[1]) for form in forms_0]
        forms_2 = [secondVariationalDerivative(form,self.params[0],self.params[2],self.params[1]) for form in forms_1]
        # Return derivatives
        return (forms_0, forms_1, forms_2)
    def add_domain(self,*forms):
        # Exit function if no forms given
        if forms == ():
            return
        # Differentiate forms
        forms_0, forms_1, forms_2 = self.differentiate_forms(*forms)
        # Add to existing forms
        self.__domain_0.extend(forms_0)
        self.__domain_1.extend(forms_1)
        self.__domain_2.extend(forms_2)
    def add_boundary(self,*forms):
        # Exit function if no forms given
        if forms == ():
            return
        # Differentiate forms
        forms_0, forms_1, forms_2 = self.differentiate_forms(*forms)
        # Add to existing forms
        self.__boundary_0.extend(forms_0)
        self.__boundary_1.extend(forms_1)
        self.__boundary_2.extend(forms_2)

class TestTrialBase: # Not intended for use except as base class
    def set_test_func(self,test_func):
        if test_func is None:
            raise TypeError
        self.__test_func = Param(test_func)
    def set_trial_func(self,trial_func=None):
        if trial_func is None:
            self.__trial_func = None
        else:
            self.__trial_func = Param(trial_func)
    @property
    def test_func(self):
        return self.__test_func
    @property
    def trial_func(self):
        return self.__trial_func

class PDE_System:
    def __init__(self,domain_PDE,boundary_PDE=None):
        self.set_domain(domain_PDE)
        self.set_boundary(boundary_PDE)
    def __repr__(self):
        if self.boundary is None:
            return f'<PDE System: {self.domain.eqn_str}>'
        else:
            return f'<PDE System: {self.domain.eqn_str}\n             {self.boundary.eqn_str}>'
    @property
    def domain(self):
        return self.__domain
    @property
    def boundary(self):
        return self.__boundary


    def set_domain(self,domain_PDE):
        if domain_PDE is None:
            raise TypeError
        if not isinstance(domain_PDE,PDE):
            raise TypeError
        self.__domain = domain_PDE
    def set_boundary(self,boundary_PDE):
        if boundary_PDE is None:
            self.__boundary = None
        elif not isinstance(boundary_PDE,PDE):
            raise TypeError
        else:
            self.__boundary = boundary_PDE

class PDE(TestTrialBase):
    def __init__(self,lhs,rhs,over='domain'):
        if not isinstance(lhs,lhsForm):
            raise TypeError
        if not isinstance(rhs,rhsForm):
            raise TypeError
        if over not in ('domain','boundary'):
            raise ValueError
        if lhs.test_func != rhs.test_func:
            print(lhs.test_func,rhs.test_func)
            raise ValueError
        self.__lhs = lhs
        self.__rhs = rhs
        self.__over = over
        self.set_trial_func(lhs.trial_func)
        self.set_test_func(lhs.test_func)
    def __repr__(self) -> str:
        return f'<PDE : {self.eqn_str}>'
    @property
    def lhs(self):
        return self.__lhs
    @property
    def rhs(self):
        return self.__rhs
    @property
    def over(self):
        return self.__over
    @property
    def eqn_str(self):
        lhs = ' + '.join([repr(form) for form in self.lhs.forms])
        rhs = ' + '.join([repr(form) for form in self.rhs.forms])
        return f'[{lhs} = {rhs}] over {self.over}'

class lhsForm(TestTrialBase):
    def __init__(self,trial_func,test_func,name=None,forms=[]):
        if not isinstance(forms,list):
            raise TypeError('\'forms\' must be a List of items of type GeneralForm.')
        self.set_trial_func(trial_func)
        self.set_test_func(test_func)
        self.name = name
        self.forms = []
        self.add_form(*forms)

    def __call__(self):
        expr = 0
        for form in self.forms:
            if not is_linear_param(form.expr,self.trial_func):
                raise ValueError(f'The form \'{form}\' of \'{self}\' is nonlinear in the trial function.')
            elif not is_linear_param(form.expr,self.test_func):
                raise ValueError(f'The form \'{form}\' of \'{self}\' is nonlinear in the test function.')
            expr += form.expr
        return uflfy(expr)

    def __repr__(self):
        return 'Untitled lhsForm' if self.name == None else f'lhsForm {self.name}'

    def add_form(self,*forms):
        for form in forms:
            if not isinstance(form,GeneralForm):
                raise TypeError('Form must be type GeneralForm.')
            self.forms.append(form)

class rhsForm(TestTrialBase):
    def __init__(self,test_func,name=None,forms=[]):
        if not isinstance(forms,list):
            raise TypeError('\'forms\' must be a List of items of type GeneralForm.')
        self.set_test_func(test_func)
        self.set_trial_func()
        self.name = name
        self.forms = []
        self.add_form(*forms)

    def __call__(self):
        expr = 0
        for form in self.forms:
            if not is_linear_param(form.expr,self.test_func):
                raise ValueError(f'The form \'{form}\' of \'{self}\' is nonlinear.')
            expr += form.expr
        return uflfy(expr)

    def __repr__(self):
        return 'Untitled rhsForm' if self.name == None else f'rhsForm {self.name}'

    def add_form(self,*forms):
        for form in forms:
            if not isinstance(form,GeneralForm):
                raise TypeError('Form must be type GeneralForm.')
            self.forms.append(form)

# Set up variable with respect to which we will take derivatives

x = [Symbol('x0'),Symbol('x1'),Symbol('x2')]

# Set up Qvector objects

nu = AbstractVector('nu')

q = QVector('q')
Dq = q.grad
Q = q.tens

p = QVector('p')
Dp = p.grad
P = p.tens

r = QVector('r')
Dr = r.grad
R = r.tens

qp = QVector('q_prev')
Dqp = qp.grad
QP = qp.tens

qpp = QVector('q_prev_prev')
Dqpp = qpp.grad
QPP = qpp.tens

qnp = QVector('q_newt_prev')
Dqnp = qnp.grad
QNP = qnp.tens

f = QVector('f')
g = QVector('g')

# END OF CODE
