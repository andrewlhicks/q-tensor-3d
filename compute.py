from sympyplus import *
from misc import getValues
from settings import const, options, userfunc

# Set up variable with respect to which we will take derivatives

x = [Symbol('x0'),Symbol('x1'),Symbol('x2')]

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

nu = AbstractVector('nu')

# Compute boundary and initial guess

def boundary():   
    n = userfunc.boundary()
    G = outerp(n,n) - (1.0/3.0) * eye(3)
    
    return uflfy(vectorfy(G))

def initialGuess():
    return uflfy(userfunc.initialGuess())

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

#  Bilinear and linear to be integrated in the domain and on the boundary

def bilinear():
    # Combine
        
    expression = (1/const.dt) * q.dot(p) + bilinearFormElastic(Dq,q,Dp,p) + (1/const.ep) * bilinearFormBulkC(Dq,q,Dp,p)

    # Convert to UFL and return
    
    return uflfy(expression)

def linear():
    # Combine
    
    expression = (1/const.dt) * qp.dot(p) + (1/const.ep) * bilinearFormBulkE(Dqp,qp,Dp,p) + linearFormForcing(Dp,p)
    
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
    return innerp(grad1,grad2)

def termL2(grad1,grad2):
    term = 0

    for ii in range(5):
        for jj in range(5):
            term += (grad1.row(ii) * E[ii]).dot(E[jj] * grad2.row(jj).T)
    
    return term

def termL3(grad1,grad2):
    term = 0

    for ii in range(5):
        for jj in range(5):
            term += (grad1.row(ii) * E[jj]).dot(E[ii] * grad2.row(jj).T)
    
    return term

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

energyElastic = Lagrangian(const.L1/2*termL1(Dq,Dq)+const.L2/2*termL2(Dq,Dq)+const.L3/2*termL3(Dq,Dq),Dq,q)
energyBulkC = Lagrangian((const.L0/2)*innerp(Q,Q),Dq,q)
energyBulkE = Lagrangian((const.L0+const.A)/2*trace(Q**2) + (const.B/3)*trace(Q**3) - (const.C/4)*trace(Q**2)**2,Dq,q)

W1 = 1
W2 = 1
S0 = (const.B + sqrt(const.B**2 + 24*const.A*const.C))/(4*const.C)
QT = Q + S0*eye(3)
Pi = eye(3) - outerp(nu,nu)
energyPlanarAnchoring = Lagrangian(W1*innerp(QT-Pi*QT*Pi,QT-Pi*QT*Pi) + W2*(innerp(QT,QT) - S0**2)**2,Dq,q)

##############

bilinearFormElastic = variationalDerivative(energyElastic,Dq,q,Dp,p)
bilinearFormBulkC = variationalDerivative(energyBulkC,Dq,q,Dp,p)
bilinearFormBulkE = variationalDerivative(energyBulkE,Dq,q,Dp,p)
bilinearFormPlanarAnchoring = variationalDerivative(energyPlanarAnchoring,Dq,q,Dp,p)

##############

def linearFormForcing(Dp,p):
    if not isinstance(Dp,AbstractVectorGradient):
        raise TypeError('Argument \'Dp\' must be type AbstractVectorGradient')
    if not isinstance(p,QVector):
        raise TypeError('Argument \'p\' must be type QVector')

    # Create a forcing if the manufactured solution is set
    
    if options.manufactured:
        n = userfunc.boundary()
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
    return const.L1*strongL1(G) + const.L2*strongL2(G) + const.L3*strongL3(G) + (1/const.ep)*(-const.A*strongA(G) - const.B*strongB(G) + const.C*strongC(G))

def strongL1(Q):
    term = zeros(3,3)
    
    for ii in range(3):
            for jj in range(3):
                for kk in range(3):
                    term[ii,jj] -= diff(Q[ii,jj],x[kk],x[kk])
    
    return term

def strongL2(Q):
    term = zeros(3,3)
        
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                term[ii,jj] -= diff(Q[ii,kk],x[jj],x[kk])
    
    return term

def strongL3(Q):  
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