from sympyplus import *
from misc import getValues
from settings import const, options, userBoundary

# Configure constants (Might want to take care of this in settings.py)

L1,L2,L3 = getValues(const,'L1 L2 L3') # Elastic constants
A,B,C = getValues(const,'A B C') # Bulk constants
L0 = getValues(const,'L0')  # Convex splitting constant
ep = getValues(const,'ep') # Epsilon
dt = getValues(const,'dt') # Time step

# Set up Qvector objects

q = QVector('q')
Dq = q.grad
Q = q.tens

qp = QVector('q_prev')
Dqp = qp.grad
QP = qp.tens

p = QVector('p')
Dp = p.grad
P = p.tens

# Compute boundary and initial guess

def boundary():
    from settings import userBoundary
    from sympy import eye
    from sympyplus import outerp, uflfy, vectorfy
    
    n = userBoundary()
    G = outerp(n,n) - (1.0/3.0) * eye(3)
    
    return uflfy(vectorfy(G))

def initialGuess():
    from settings import userInitialGuess
    from sympyplus import uflfy
    
    return uflfy(userInitialGuess())

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

#  Bilinear and linear to be integrated in the domain and on the boundary

def bilinear():
    from sympyplus import uflfy
    
    # Combine
        
    expression = (1/dt) * q.dot(p) + bilinearFormElastic(Dq,q,Dp,p) + (1/ep) * bilinearFormBulkC(Dq,q,Dp,p)

    # Convert to UFL and return
    
    return uflfy(expression)

def linear():
    from sympyplus import uflfy
    
    # Combine
    
    expression = (1/dt) * qp.dot(p) + (1/ep) * bilinearFormBulkE(Dqp,qp,Dp,p) + linearFormForcing(Dp,p)
    
    # Convert to UFL and return
    
    return uflfy(expression)

def bilinearOnBoundary():
    pass

def linearOnBoundary():
    pass

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

def termL1(grad1,grad2):
    from sympyplus import innerp
    
    return innerp(grad1,grad2)

def termL2(grad1,grad2):
    from sympyplus import E
    
    term = 0

    for ii in range(5):
        for jj in range(5):
            term += (grad1.row(ii) * E[ii]).dot(E[jj] * grad2.row(jj).T)
    
    return term

def termL3(grad1,grad2):
    from sympyplus import E
    
    term = 0

    for ii in range(5):
        for jj in range(5):
            term += (grad1.row(ii) * E[jj]).dot(E[ii] * grad2.row(jj).T)
    
    return term

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

energyElastic = Lagrangian(L1/2*termL1(Dq,Dq)+L2/2*termL2(Dq,Dq)+L3/2*termL3(Dq,Dq),Dq,q)
energyBulkC = Lagrangian((L0/2)*innerp(Q,Q),Dq,q)
energyBulkE = Lagrangian((L0+A)/2*trace(Q**2) + (B/3)*trace(Q**3) - (C/4)*trace(Q**2)**2,Dq,q)

##############

bilinearFormElastic = variationalDerivative(energyElastic,Dq,q,Dp,p)
bilinearFormBulkC = variationalDerivative(energyBulkC,Dq,q,Dp,p)
bilinearFormBulkE = variationalDerivative(energyBulkE,Dq,q,Dp,p)

##############

def linearFormForcing(Dp,p):
    if not isinstance(Dp,AbstractVectorGradient):
        raise TypeError('Argument \'Dp\' must be type AbstractVectorGradient')
    if not isinstance(p,QVector):
        raise TypeError('Argument \'p\' must be type QVector')

    # Create a forcing if the manufactured solution is set
    
    if options['manufactured']:
        n = userBoundary()
        G = outerp(n,n) - (1.0/3.0) * eye(3)
        F = strongForm(G)
        f = vectorfy(F)
    else:
        f = Matrix([0,0,0,0,0])
    
    # Combine and return
    
    return f.dot(p)

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

def strongForm(G): # plugs G into the strong form PDE
    from settings import const

    A, B, C, ep, L1, L2, L3 = getValues(const,'A B C ep L1 L2 L3')
    
    return L1*strongL1(G) + L2*strongL2(G) + L3*strongL3(G) + (1/ep)*(-A*strongA(G) - B*strongB(G) + C*strongC(G))

def strongL1(Q):
    from sympy import symbols, zeros, diff
    x0,x1,x2 = symbols('x0 x1 x2')
    x = [x0,x1,x2]
    
    term = zeros(3,3)
    
    for ii in range(3):
            for jj in range(3):
                for kk in range(3):
                    term[ii,jj] -= diff(Q[ii,jj],x[kk],x[kk])
    
    return term

def strongL2(Q):
    from sympy import symbols, zeros, diff
    x0,x1,x2 = symbols('x0 x1 x2')
    x = [x0,x1,x2]
    
    term = zeros(3,3)
        
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                term[ii,jj] -= diff(Q[ii,kk],x[jj],x[kk])
    
    return term

def strongL3(Q):
    from sympy import symbols, zeros, diff
    x0,x1,x2 = symbols('x0 x1 x2')
    x = [x0,x1,x2]
    
    term = zeros(3,3)
        
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                term[ii,jj] -= diff(Q[ii,kk],x[kk],x[jj])
    
    return term

def strongA(Q):
    return Q

def strongB(Q):
    return Q*Q

def strongC(Q):
    from sympyplus import innerp
    
    return innerp(Q,Q)*Q

# END OF CODE